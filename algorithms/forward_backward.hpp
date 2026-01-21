#pragma once

#include <vector>
#include <array>
#include <cmath>

#include "../utils/structs_consts_functions.hpp"

using namespace std;


/**
 * @brief Forward algoritam sa skaliranjem kako bi izbjegli padanje
 * vrijednosti alpha u 0.
 * 
 * @param O Sekvenca opažanja
 * @param hmm HMM parametri
 * @param alpha Matrica za pohranu forward varijabli
 * @param c Vektor skalirajućih faktora
 * 
 * @return double Log-vjerojatnost sekvence
 */
double forward_scaled(
    const vector<int>& O,
    const HMM& hmm,
    vector<array<double, NSTATE>>& alpha,
    vector<double>& c
);


/**
 * @brief Forward algoritam sa skaliranjem uz masku dozvoljenih stanja po vremenu.
 *
 * @param O Sekvenca opažanja
 * @param hmm HMM parametri
 * @param state_mask Maska dozvoljenih stanja (0/1) za svaki t
 * @param alpha Matrica za pohranu forward varijabli
 * @param c Vektor skalirajućih faktora
 *
 * @return double Log-vjerojatnost sekvence
 */
double forward_scaled_masked(
    const vector<int>& O,
    const HMM& hmm,
    const vector<array<double, NSTATE>>& state_mask,
    vector<array<double, NSTATE>>& alpha,
    vector<double>& c
);


/**
 * @brief Backward algoritam sa skaliranjem kako bi izbjegli padanje
 * vrijednosti beta u 0. Koristi skalirajuće faktore iz forward algoritma.
 * 
 * @param O Niz opažanja
 * @param hmm HMM parametri
 * @param c Vektor skalirajućih faktora iz forward algoritma
 * @param beta Matrica za pohranu backward varijabli
 */
void backward_scaled(
    const vector<int>& O,
    const HMM& hmm,
    const vector<double>& c,
    vector<array<double, NSTATE>>& beta
);

/**
 * @brief Backward algoritam sa skaliranjem uz masku dozvoljenih stanja po vremenu.
 *
 * @param O Niz opažanja
 * @param hmm HMM parametri
 * @param state_mask Maska dozvoljenih stanja (0/1) za svaki t
 * @param c Vektor skalirajućih faktora iz forward algoritma
 * @param beta Matrica za pohranu backward varijabli
 */
void backward_scaled_masked(
    const vector<int>& O,
    const HMM& hmm,
    const vector<array<double, NSTATE>>& state_mask,
    const vector<double>& c,
    vector<array<double, NSTATE>>& beta
);
