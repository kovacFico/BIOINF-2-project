#include "./decode.hpp"
#include "../postprocesing/decoded_postprocesing.hpp"


vector<double> compute_posterior_c(const vector<int>& Oseg, const HMM& hmm) {
    vector<array<double, NSTATE>> alpha, beta;
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


vector<int> decode_hysteresis(const vector<double>& posterior_c, double enter_th, double exit_th) {
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


vector<CpgRegion> process_window(
    const vector<int>& O,
    const string& sequence,
    const HMM& hmm,
    int start_d,
    int end_d,
    int T,
    double POST_ENTER, 
    double POST_EXIT, 
    double POST_TRIM, 
    int OVERLAP
) {
    vector<int> Oseg(O.begin() + start_d, O.begin() + end_d);

    auto posterior = compute_posterior_c(Oseg, hmm);
    auto states    = decode_hysteresis(posterior, POST_ENTER, POST_EXIT);

    vector<CpgRegion> islands;
    extract_cpg_islands(islands, states);

    // Globalni pomak za početak segmenta
    int base_shift = start_d;
    for (auto& r : islands) {
        r.start += base_shift;
        r.end   += base_shift;
    }

    int keep_left_d  = (start_d == 0) ? start_d : start_d + OVERLAP / 2;
    int keep_right_d = (end_d == T)   ? end_d  : end_d  - OVERLAP / 2;

    keep_and_clip(islands, keep_left_d + 1, keep_right_d + 1); // +1 zbog 1-based koordinata
    trim_islands_with_posterior(islands, posterior, base_shift, POST_TRIM);
    filter_lenght_and_merge_close_islands(islands);
    filter_by_content(sequence, islands);

    cout << "Window " << start_d << "-" << end_d
         << " (bp keep " << keep_left_d + 1 << "-" << keep_right_d + 1
         << "), islands=" << islands.size() << "\n";

    return islands;
}
