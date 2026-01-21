#include "../hmm/hmm_io.hpp"
#include "../hmm/hmm.hpp"
#include "../utils/structs_consts_functions.hpp"
#include "../train_functions/train_func.hpp"
#include "../algorithms/baum_welch.hpp"


int main() {
    HMM hmm;

    if (ifstream("../output/trained_hmm_params.txt")) {
        hmm = load_hmm("../output/trained_hmm_params.txt");
    } else {
        hmm = load_hmm("../output/init_hmm_params.txt");
    }

    string s;
    vector<lowerCaseRegions> lc;
    get_chromosome_and_lowercase_regions("../output/" + to_string(hmm.chromosome) + "_train_chr.txt", s, lc);
    
    cout << "Učitana sekvenca za kromosom " << hmm.chromosome << " dužine " << s.size() << endl;

    vector<CpgRegion> coords_chr_orig = load_all_or_selected_coords(hmm.chromosome);
    vector<CpgRegion> coords_chr_comp = map_orig_coords_to_compressed(s, lc, coords_chr_orig);

    // ------- SEMI-SUPERVIZIJA: maska dozvoljenih stanja -------
    vector<vector<int>> sequences;
    vector<vector<array<double, NSTATE>>> masks;
    build_masked_sequences(s, coords_chr_comp, sequences, masks);
    
    cout << "Izgrađene " << sequences.size() << " trening sekvence sa maskama.\n";

    // ------- Baum-Welch na mini-sekvencama -------
    double prev_ll = -1e100;
    for (int iter = 0; iter < 10; iter++) {
        double ll = 0.0;
        baum_welch_iteration_multi_masked(sequences, masks, hmm, ll);
        
        cout << "Iter " << iter << " logL = " << ll << endl;

        if (fabs(ll - prev_ll) < 1e-3) break;
        prev_ll = ll;
    }

    save_hmm(hmm, "../output/trained_hmm_params.txt");

    return 0;
}
