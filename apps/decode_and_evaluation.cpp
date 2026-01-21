#include "../algorithms/forward_backward.hpp"
#include "../hmm/hmm_io.hpp"
#include "../postprocesing/decoded_postprocesing.hpp"
#include "../evaluation/evaluation.hpp"
#include "../utils/structs_consts_functions.hpp"

// Zadrži samo predikcije koje se preklapaju s intervalom zadržavanja u baznim koordinatama
// i skrati ih na taj interval (da ne uzimamo rubove prozora).
static void keep_and_clip(vector<CpgRegion>& islands, int keep_start_bp, int keep_end_bp) {
    vector<CpgRegion> out;
    out.reserve(islands.size());

    for (auto r : islands) {
        int s = max(r.start, keep_start_bp);
        int e = min(r.end, keep_end_bp);
        if (e >= s) out.push_back({s, e, 0});
    }

    islands.swap(out);
}

// Izračunaj posteriornu vjerojatnost P(Z_t = CpG | Oseg) za svaki dinukleotid t koristeći forward/backward.
static vector<double> compute_posterior_c(const vector<int>& Oseg, const HMM& hmm) {
    vector<array<double, NSTATE>> alpha;
    vector<array<double, NSTATE>> beta;
    vector<double> c;

    forward_scaled(Oseg, hmm, alpha, c);
    backward_scaled(Oseg, hmm, c, beta);

    vector<double> posterior(Oseg.size(), 0.0);
    for (size_t t = 0; t < Oseg.size(); t++) {
        double norm = 0.0;
        for (int i = 0; i < NSTATE; i++) {
            norm += alpha[t][i] * beta[t][i];
        }
        if (norm > 0.0) {
            posterior[t] = (alpha[t][1] * beta[t][1]) / norm;
        }
    }

    return posterior;
}

// Pretvori posterior u tvrda stanja uz histerezu (prag ulaska/izlaska) radi stabilnih CpG otoka bez treperenja.
static vector<int> decode_hysteresis(const vector<double>& posterior_c, double enter_th, double exit_th) {
    vector<int> states(posterior_c.size(), 0);
    bool in_cpg = false;

    for (size_t t = 0; t < posterior_c.size(); t++) {
        if (!in_cpg && posterior_c[t] >= enter_th) {
            in_cpg = true;
        } else if (in_cpg && posterior_c[t] < exit_th) {
            in_cpg = false;
        }
        states[t] = in_cpg ? 1 : 0;
    }

    return states;
}

// Suzi granice CpG otoka uklanjanjem rubova gdje je posterior ispod praga trim_th.
static void trim_islands_with_posterior(
    vector<CpgRegion>& islands,
    const vector<double>& posterior_c,
    int base_shift,
    double trim_th
) {
    vector<CpgRegion> trimmed;
    trimmed.reserve(islands.size());

    const int max_d = static_cast<int>(posterior_c.size());

    for (const auto& r : islands) {
        int local_start_bp = r.start - base_shift;
        int local_end_bp = r.end - base_shift;

        int d_start = local_start_bp;
        int d_end = local_end_bp - 1;

        if (d_start < 1) d_start = 1;
        if (d_end > max_d) d_end = max_d;

        while (d_start <= d_end && posterior_c[(size_t)(d_start - 1)] < trim_th) {
            d_start++;
        }
        while (d_end >= d_start && posterior_c[(size_t)(d_end - 1)] < trim_th) {
            d_end--;
        }

        if (d_start <= d_end) {
            int new_start_bp = d_start + base_shift;
            int new_end_bp = d_end + 1 + base_shift;
            if (new_end_bp >= new_start_bp) {
                trimmed.push_back({new_start_bp, new_end_bp, r.chromosome});
            }
        }
    }

    islands.swap(trimmed);
}


void load_chr_seq_to_dinuc_vector(vector<int>& O, int chr_number) {
    ifstream in("../output/" + to_string(chr_number) + "_test_chr.txt");
    if (!in) {
        cerr << "Ne mogu otvoriti test fajl za kromosom " << chr_number << "\n";
        exit(1);
    }
    string s;
    getline(in, s);

    O.reserve(s.size() - 1);
    for (size_t i = 1; i < s.size(); i++) {
        int x = di_index(s[i - 1], s[i]);
        if (x != -1) O.push_back(x);
    }
}

int main() {
    HMM hmm = load_hmm("../output/trained_hmm_params.txt");

    if (hmm.chromosome < 17) hmm.chromosome = 17;

    vector<int> O;
    ifstream in("../output/" + to_string(hmm.chromosome) + "_test_chr.txt");
    if (!in) {
        cerr << "Ne mogu otvoriti test fajl za kromosom " << hmm.chromosome << "\n";
        exit(1);
    }
    string s;
    getline(in, s);

    O.reserve(s.size() - 1);
    for (size_t i = 1; i < s.size(); i++) {
        int x = di_index(s[i - 1], s[i]);
        if (x != -1) O.push_back(x);
    }

    cout << "Učitana sekvenca za kromosom " << hmm.chromosome
              << " (baze=" << s.size()
              << ", dinukleotidi=" << O.size() << ")\n";

    
    // Windowing parametri (dinukleotidi)
    const int WIN = 5'000'000;       // broj dinukleotida po prozoru
    const int OVERLAP = 50'000;      // preklapanje (dinukleotidi)
    const int STEP = WIN - OVERLAP;  // pomak prozora
    const double POST_ENTER = 0.60;
    const double POST_EXIT = 0.40;
    const double POST_TRIM = 0.42;

    vector<CpgRegion> predicted_all;
    predicted_all.reserve(20000);

    int T = (int)O.size();
    for (int start_d = 0; start_d < T; start_d += STEP) {
        int end_d = min(start_d + WIN, T);
        int len_d = end_d - start_d;
        if (len_d < 2) break;

        // Segment opažanja
        vector<int> Oseg;
        Oseg.assign(O.begin() + start_d, O.begin() + end_d);

        // Posterior + dekodiranje s histerezom na segmentu.
        vector<double> posterior_c = compute_posterior_c(Oseg, hmm);
        vector<int> states = decode_hysteresis(posterior_c, POST_ENTER, POST_EXIT);

        // CpG otoci u baznim koordinatama segmenta (zbog mapiranja end+1 u extract_cpg_islands).
        vector<CpgRegion> islands_seg;
        extract_cpg_islands(islands_seg, states);

        // Globalni offset:
        // start_d je 0-based dinukleotid index u globalnom O.
        // Dinukleotid t (0-based) pokriva baze (t+1, t+2) u globalnim bazama.
        // Naš extract vraća bazne koordinate 1-based unutar segmenta.
        // Zato global shift u bazama = start_d
        // (jer segment baza start = start_d+1)
        int base_shift = start_d;

        for (auto& r : islands_seg) {
            r.start += base_shift;
            r.end   += base_shift;
        }

        // Zadrži samo “sredinu” prozora da rubovi ne rade artefakte:
        // keep interval u bazama:
        // - lijevi rub: za prvi prozor zadržimo od početka
        // - inače odrežemo prvih OVERLAP/2 dinukleotida (~ baza)
        // - desni rub: za zadnji prozor zadržimo do kraja
        int keep_left_d  = (start_d == 0) ? start_d : start_d + OVERLAP / 2;
        int keep_right_d = (end_d == T)   ? end_d  : end_d  - OVERLAP / 2;

        // Dinukleotidi -> baze: [d_left, d_right] dinukleotidi pokrivaju baze [d_left+1, d_right+1]
        int keep_start_bp = keep_left_d + 1;
        int keep_end_bp   = keep_right_d + 1;

        keep_and_clip(islands_seg, keep_start_bp, keep_end_bp);
        trim_islands_with_posterior(islands_seg, posterior_c, base_shift, POST_TRIM);
        filter_lenght_and_merge_close_islands(islands_seg);
        filter_by_content(s, islands_seg);

        // dodaj u globalnu listu
        predicted_all.insert(predicted_all.end(), islands_seg.begin(), islands_seg.end());

        cerr << "Window " << start_d << "-" << end_d
                  << " (bp keep " << keep_start_bp << "-" << keep_end_bp
                  << "), islands=" << islands_seg.size() << "\n";
    }


    move_predicted_based_on_lowercase(predicted_all, hmm.chromosome);

    vector<CpgRegion> true_islands;
    load_true_islands(true_islands, hmm.chromosome);
    island_based_evaluation(predicted_all, true_islands);
    base_pair_evaluation(predicted_all, true_islands);

    save_hmm(hmm, "../output/trained_hmm_params.txt");

    return 0;
}
