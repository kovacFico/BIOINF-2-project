#include "../algorithms/decode.hpp"
#include "../hmm/hmm_io.hpp"
#include "../hmm/hmm.hpp"
#include "../postprocesing/decoded_postprocesing.hpp"
#include "../evaluation/evaluation.hpp"
#include "../utils/structs_consts_functions.hpp"


// Windowing parametri (dinukleotidi)
const int WINDOW = 5'000'000;       // broj dinukleotida po prozoru
const int OVERLAP = 50'000;      
const int STEP = WINDOW - OVERLAP; 
const double POST_ENTER = 0.60;
const double POST_EXIT = 0.40;
const double POST_TRIM = 0.42;

/* 
 * @brief Predikcija CpG otoka pomoću treniranog HMM-a uz prozorsku obradu
 *        cijelog kromosoma i evaluaciju rezultata.
 *
 * Program učitava trenirani skriveni Markovljev model (HMM) i sekvencu
 * odabranog kromosoma, pretvara baznu sekvencu u dinukleotidna opažanja
 * te provodi dekodiranje u preklapajućim prozorima.
 *
 * Za svaki prozor:
 *  - računa se posterior CpG stanja (forward-backward)
 *  - provodi se dekodiranje s histerezom
 *  - ekstrahiraju se CpG otoci
 *  - uklanjaju se rubni artefakti i provodi sadržajno filtriranje
 *
 * Predviđeni CpG otoci se zatim:
 *  - korigiraju prema lowercase regijama
 *  - uspoređuju s referentnim anotacijama
 *  - evaluiraju na razini otoka i baznih parova
 *
 * Parametri prozora (veličina, preklapanje) i pragovi posteriora
 * definirani su kao globalne konstante.
 *
 * @note Koordinate CpG otoka su izražene u 1-based baznim koordinatama.
 * @note Dinukleotidni indeksi su 0-based.
 */
int main() {
    HMM hmm = load_hmm("../output/trained_hmm_params.txt");
    if (hmm.chromosome < 17) hmm.chromosome = 17;

    vector<int> O;
    string s;
    load_chr_seq_to_dinuc_vector(O, s, hmm.chromosome);

    cout << "Učitana sekvenca za kromosom " << hmm.chromosome
         << " (baze=" << s.size()
         << ", dinukleotidi=" << O.size() << ")\n";

    vector<CpgRegion> predicted_all;
    predicted_all.reserve(20000);

    int T = (int)O.size();
    for (int start_d = 0; start_d < T; start_d += STEP) {
        int end_d = min(start_d + WINDOW, T);
        int len_d = end_d - start_d;
        if (len_d < 2) break;

        auto islands = process_window(
            O, s, hmm, start_d, end_d, T,
            POST_ENTER, POST_EXIT, POST_TRIM, OVERLAP
        );

        predicted_all.insert(
            predicted_all.end(),
            islands.begin(),
            islands.end()
        );
    }

    move_predicted_based_on_lowercase(predicted_all, hmm.chromosome);

    vector<CpgRegion> true_islands = load_all_or_selected_coords(hmm.chromosome);
    
    island_based_evaluation(predicted_all, true_islands);
    base_pair_evaluation(predicted_all, true_islands);

    save_hmm(hmm, "../output/trained_hmm_params.txt");

    return 0;
}
