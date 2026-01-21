#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

#include "../utils/structs_consts_functions.hpp"

using namespace std;


/**
 * @brief Evaluira predviđene CpG otoke u odnosu na stvarne otoke.
 * Island based evaluacija: broji se koliko je predviđenih otoka točno
 * predviđeno (TP), koliko je netočno predviđeno (FP) i koliko
 * stvarnih otoka nije predviđeno (FN). Ne računa se boundary precision
 * 
 * @param predicted Vektor predviđenih CpG otoka
 * @param truth Vektor stvarnih CpG otoka
 * 
 * Napomena: moguća prenamjena na base pair evaluaciju ako bude potrebno
 */
void island_based_evaluation(const vector<CpgRegion>& predicted, const vector<CpgRegion>& truth);


/**
 * @brief Evaluira predviđene CpG otoke na razini baza u odnosu na stvarne otoke.
 * Base-pair based evaluacija: računa se broj točno predviđenih baza (TP),
 * netočno predviđenih baza (FP) i stvarnih baza koje nisu predviđene (FN).
 * 
 * @param predicted Vektor predviđenih CpG otoka
 * @param truth Vektor stvarnih CpG otoka
 */
void base_pair_evaluation(const vector<CpgRegion>& predicted, const vector<CpgRegion>& truth);