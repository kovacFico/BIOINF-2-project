#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include "../utils/structs_consts_functions.hpp"

using namespace std;



// Minimalna dužina CpG otoka kako bi spriječili šum, tj
// kratke predviđene otoke ignoriramo
constexpr int MIN_CPG_LEN = 250;
constexpr int MERGE_DISTANCE = 50;
constexpr double MIN_GC_CONTENT = 0.44;
constexpr double MIN_CPG_OE = 0.50;

/**
 * @brief Izvlači CpG otoke iz binarne sekvence stanja
 * pod pretpostavkom da je stanje 1 CpG otok, a stanje 0 ne-CpG regija. Također
 * primjenjuje minimalnu dužinu CpG otoka kako bi se spriječio šum, tj. kratke
 * predviđene otoke ignorira.
 * 
 * @param islands Referenca na vektor za pohranu CpG otoka
 * @param states Sekvenca stanja (0=B, 1=C)
 */
void extract_cpg_islands(vector<CpgRegion>& islands, vector<int>& states);


/**
 * @brief Pomiče predviđene CpG otoke na temelju malih slova u genomu
 * 
 * @param predicted Vektor predviđenih CpG otoka
 * @param chr_number Broj kromosoma
 */
void move_predicted_based_on_lowercase(vector<CpgRegion>& predicted, int chr_number);


/**
 * @brief Filtrira CpG otoke prema MIN_CPG_LEN dužini i spaja otoke koji su
 * udaljeni najviše MERGE_DISTANCE baza.
 * 
 * @param islands Referenca na vektor CpG otoka
 */
void filter_lenght_and_merge_close_islands(vector<CpgRegion>& islands);


/**
 * @brief Filtrira CpG otoke prema GC sadržaju i CpG O/E omjeru.
 *
 * @param sequence Referenca na sekvencu kromosoma (1-based logika u CpgRegion)
 * @param islands Referenca na vektor CpG otoka
 */
void filter_by_content(const string& sequence, vector<CpgRegion>& islands);



/**
 * @brief Učitava točne CpG otoke iz coords datoteke. 
 * 
 * @param islands Referenca na vektor za pohranu CpG otoka
 * @param chr_number Broj kromosoma čije se otoke učitava
 * 
 */
void load_true_islands(vector<CpgRegion>& islands, int chr_number);

