#include "./hmm_io.hpp"
#include "./hmm.hpp"

#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

HMM load_hmm(const string &filename) {
    ifstream in(filename);
    if (!in) {
        cerr << "Ne mogu otvoriti HMM parametre\n";
        exit(1);
    }

    HMM hmm;
    string tmp;

    if (filename == "../output/trained_hmm_params.txt") {
        in >> tmp >> hmm.chromosome;  
    } else {
        hmm.chromosome = 1;
    }

    in >> tmp >> tmp;                
    for (int k = 0; k < NSYM; k++)
        in >> hmm.B[0][k];

    in >> tmp >> tmp;                 
    for (int k = 0; k < NSYM; k++)
        in >> hmm.B[1][k];

    in >> tmp;                       
    in >> tmp >> hmm.A[0][0];         
    in >> tmp >> hmm.A[0][1];         
    in >> tmp >> hmm.A[1][1];         
    in >> tmp >> hmm.A[1][0];         

    in >> tmp >> hmm.pi[0] >> hmm.pi[1];

    return hmm;
}


void save_hmm(const HMM& hmm, const string& filename) {
    ofstream out(filename);
    out << fixed << setprecision(8);

    if (filename == "../output/trained_hmm_params.txt") {
        out << "Kromosom: " << hmm.chromosome + 1 << "\n";
    }

    out << "Emisije B: ";
    for (int k = 0; k < NSYM; k++) out << hmm.B[0][k] << " ";
    out << "\n";

    out << "Emisije C: ";
    for (int k = 0; k < NSYM; k++) out << hmm.B[1][k] << " ";
    out << "\n";

    out << "Tranzicije:\n";
    out << "B->B " << hmm.A[0][0] << "\n";
    out << "B->C " << hmm.A[0][1] << "\n";
    out << "C->C " << hmm.A[1][1] << "\n";
    out << "C->B " << hmm.A[1][0] << "\n";

    if (filename == "../output/trained_hmm_params.txt") {
        out << "PI: " << hmm.pi[0] << " " << hmm.pi[1] << "\n";
    } else {
        out << "PI: 0.9 0.1\n";
    }
}
