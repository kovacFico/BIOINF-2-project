#include "./train_func.hpp"


void get_chromosome_and_lowercase_regions(const string& filename, string& seq, vector<lowerCaseRegions>& lc) {
    ifstream in(filename);
    if (!in) {
        cerr << "Ne mogu otvoriti: " << filename << endl;
        exit(1);
    }

    getline(in, seq);
    int s, e;
    while (in >> s >> e) {
        lc.push_back({s, e});
    }
}


vector<int> seq_to_dinuc(const string& s, int start, int end) {
    vector<int> O;
    int L = end - start + 1;

    if (L < 2) return O;

    O.reserve(L - 1);
    for (int pos = start + 1; pos <= end; pos++) {
        char prev = s[(size_t)(pos - 2)]; 
        char cur  = s[(size_t)(pos - 1)]; 
        int x = di_index(prev, cur);
        if (x != -1) O.push_back(x);
    }

    return O;
}


int map_to_comp(const vector<lowerCaseRegions>& lc, int orig_pos) {
    long long removed = 0;

    for (const auto& r : lc) {
        if (r.start >= orig_pos) break;
        int end = min(r.end, orig_pos - 1);
        if (end >= r.start) {
            removed += (end - r.start + 1);
        }
    }

    return orig_pos - removed;
}


vector<CpgRegion> map_orig_coords_to_compressed(
    const string& seq,
    const vector<lowerCaseRegions>& lc,
    const vector<CpgRegion>& orig_coords
) {
    vector<CpgRegion> comp_coords;
    comp_coords.reserve(orig_coords.size());

    for (const auto& r : orig_coords) {
        int c_start = map_to_comp(lc, r.start);
        int c_end = map_to_comp(lc, r.end);
        if (c_start < 1) c_start = 1;
        if (c_end > (int)seq.size()) c_end = (int)seq.size();
        if (c_start <= c_end) {
            comp_coords.push_back({c_start, c_end, r.chromosome});
        }
    }
        
    return comp_coords;
}


void build_masked_sequences(
    const string& s,
    const vector<CpgRegion>& coords_chr,
    vector<vector<int>>& sequences,
    vector<vector<array<double, NSTATE>>>& masks
) {
    const int NEG_MARGIN = 200;
    const int CHUNK_D = 1'000'000;

    /*
     * 1. Priprema pomoćnih nizova po bazama
     *
     * base_is_cpg[pos]   = 1 ako je baza unutar poznate CpG regije
     * base_near_cpg[pos] = 1 ako je baza unutar ± NEG_MARGIN od CpG regije
     *
     * Indeksiranje je logički 1-based kako bi se podudaralo s biološkim
     * koordinatama i postojećim kodom.
     */
    vector<char> base_is_cpg(s.size() + 1, 0);
    vector<char> base_near_cpg(s.size() + 1, 0);

    for (const auto& r : coords_chr) {
        int a = max(1, r.start);
        int b = min((int)s.size(), r.end);
        for (int pos = a; pos <= b; pos++) {
            base_is_cpg[pos] = 1;
        }

        int na = max(1, r.start - NEG_MARGIN);
        int nb = min((int)s.size(), r.end + NEG_MARGIN);
        for (int pos = na; pos <= nb; pos++) {
            base_near_cpg[pos] = 1;
        }
    }

    /*
     * 2. Chunkiranje sekvence
     *
     * HMM opažanja su dinukleotidi, pa je ukupan broj opažanja: T_full = |s| - 1
    */
    int T_full = int(s.size()) - 1;

    for (int start_d = 0; start_d < T_full; start_d += CHUNK_D) {
        int end_d = min(start_d + CHUNK_D, T_full);

        // odrđivanje dinukleotida u chunku
        int start_bp = start_d + 1;
        int end_bp = end_d + 1;
        auto O = seq_to_dinuc(s, start_bp, end_bp);
        if (O.size() < 2) continue;

        vector<array<double, NSTATE>> mask;
        mask.reserve(O.size());


        /* 3. Izgradnja maske dozvoljenih stanja po dinukleotidu */
        for (int d = start_d; d < end_d; d++) {
            int b1 = d + 1;
            int b2 = d + 2;
            if (b2 > (int)s.size()) break;

            // Ako je barem jedna baza unutar CpG regije -> CpG stanje
            if (base_is_cpg[b1] || base_is_cpg[b2]) {
                mask.push_back({0.0, 1.0});
            
            // Ako je dinukleotid daleko od svih CpG regija -> non-CpG
            } else if (!base_near_cpg[b1] && !base_near_cpg[b2]) {
                mask.push_back({1.0, 0.0});

            // Prijelazna zona
            } else {
                mask.push_back({1.0, 1.0});
            }
        }

        // provjera ima li svaki dinukleoid jednu masku
        if (mask.size() == O.size()) {
            sequences.push_back(move(O));
            masks.push_back(move(mask));
        }
    }
}
