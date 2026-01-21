#pragma once

#include "../utils/structs_consts_functions.hpp"
#include "../algorithms/forward_backward.hpp"

using namespace std;


/**
 * @brief Izvršava jednu iteraciju Baum-Welch algoritma (EM) za skup sekvenci
 * uz korištenje maski dozvoljenih stanja po t.
 *
 * @param sequences Vektor sekvenci opažanja (dinukleotidi)
 * @param state_masks Vektor maski dozvoljenih stanja (0/1) za svaku sekvencu
 * @param hmm HMM model čiji se parametri ažuriraju
 * @param ll Referenca na log-vjerojatnost koja se ažurira
 *
 * @return double Ažurirana log-vjerojatnost svih sekvenci
 */
double baum_welch_iteration_multi_masked(
    const vector<vector<int>>& sequences,
    const vector<vector<array<double, NSTATE>>>& state_masks,
    HMM& hmm,
    double& ll
);
