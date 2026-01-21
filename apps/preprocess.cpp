#include "../preprocesing/genome_preprocesing.hpp"


/**
 * @brief Predobrada genomskih podataka i priprema ulaznih datoteka za HMM.
 *
 * Ova funkcija provodi cjelokupni postupak predobrade:
 * - Učitava poznate koordinate pozitivnih CpG otoka.
 * - Učitava i čisti pozadinsku (background) genomsku sekvencu.
 * - Dijeli genom na pojedinačne kromosome.
 * - Bilježi lowercase (maskirane) regije za svaki kromosom.
 * - Sprema per-kromosomske sekvence i pripadne metapodatke.
 *
 * Izlazne datoteke:
 *  - <chr>_train_chr.txt / <chr>_test_chr.txt :
 *      prvi red = sekvenca kromosoma (samo uppercase),
 *      sljedeći redovi = intervali lowercase regija (1-based)
 *  - clean_positive.txt   : sekvence pozitivnih CpG otoka
 *  - clean_background.txt : pozadinski genom bez CpG regija
 *  - coords.txt           : koordinate CpG otoka (chr, start, end)
 *
 * Ova aplikacija se pokreće jednom prije inicijalizacije i treniranja HMM-a.
 *
 * @note Putanje do datoteka su trenutno zadane u kodu (hardcoded).
 */
int main() {
    const string output_dir = "../output";
    const string genome_path = "../data/ncbi_dataset/ncbi_dataset/data/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna";
    const int NUM_CHROMOSOMES = 22;
    
    string background;
    vector<CpgRegion> coords;
    vector<string> positive_cpg = load_positive_cpg("../data/test.txt", coords);
    background = load_background(genome_path, coords);

    vector<ofstream> chromosome_out_files;
    ofstream out1, out2, coords_out;
    open_output_files(NUM_CHROMOSOMES, output_dir, chromosome_out_files, out1, out2, coords_out);

    long chromosome_length = 0;
    for (int chr = 1; chr <= NUM_CHROMOSOMES; chr++) {
        vector<lowerCaseRegions> lowercaseCoords;
        string chromosome = load_chromosome(genome_path, lowercaseCoords, chr);

        cout << "Duzina kromosoma " << chr << ": " << chromosome.size() << endl;

        chromosome_length += chromosome.size();
        chromosome_out_files[chr - 1] << chromosome << "\n";

        for (const auto &l : lowercaseCoords) 
            chromosome_out_files[chr - 1] << l.start << " " << l.end << "\n";
    }

    cout << "Broj pozitivnih CpG otoka: " << positive_cpg.size() << endl;
    cout << "Duzina originalnog kromosoma: " << chromosome_length << endl;
    cout << "Duzina backgrounda nakon ciscenja: " << background.size() << endl;

    for (const auto &s : positive_cpg) out1 << s << "\n";
    for (const auto &c : coords) coords_out << c.chromosome << " " << c.start << " " << c.end << "\n";
    out2 << background;


    for (int chr = 0; chr < NUM_CHROMOSOMES; chr++) chromosome_out_files[chr].close();
    out1.close();
    out2.close();
    coords_out.close();

    return 0;
}