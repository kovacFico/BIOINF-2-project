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
 * @brief Zadržava samo CpG otoke koji se preklapaju s danim intervalom baza
 * te ih skraćuje (clipa) na granice tog intervala.
 *
 * Ova funkcija se koristi za uklanjanje rubnih artefakata kod obrade
 * preklapajućih prozora – zadržava se samo "pouzdana" sredina prozora, kako
 * bi se izbjegli efekti rubova, te prekidanje otoka na granicama prozora.
 *
 * @param islands Vektor CpG otoka u baznim koordinatama (modificira se in-place)
 * @param keep_start_bp Početna bazna koordinata intervala zadržavanja (uključivo)
 * @param keep_end_bp Završna bazna koordinata intervala zadržavanja (uključivo)
 */
void keep_and_clip(vector<CpgRegion>& islands, int keep_start_bp, int keep_end_bp);


/**
 * @brief Sužava granice CpG otoka uklanjanjem rubnih baza gdje je
 * posteriorna vjerojatnost CpG stanja ispod zadanog praga.
 *
 * Trim se provodi u dinukleotidnom prostoru, ali se rezultati vraćaju
 * u baznim koordinatama. Otoci koji nakon trimanja postanu prazni
 * se uklanjaju.
 *
 * @param islands Vektor CpG otoka u baznim koordinatama (modificira se in-place)
 * @param posterior Vektor posteriornih vjerojatnosti CpG stanja po dinukleotidu
 * @param base_shift Globalni pomak baza (početak segmenta u baznim koordinatama)
 * @param trim_th Prag posteriora ispod kojeg se rubovi otoka uklanjaju
 */
void trim_islands_with_posterior(
    vector<CpgRegion>& islands,
    const vector<double>& posterior_c,
    int base_shift,
    double trim_th
);


/**
 * @brief Učitava sekvencu zadanog kromosoma iz datoteke i pretvara je u
 * vektor dinukleotidnih opažanja. Sekvenca se učitava kao string baza (A,C,G,T), 
 * dok se dinukleotidi kodiraju kao cjelobrojni indeksi pomoću funkcije di_index().
 *
 * @param dinucs Referenca na vektor u koji se spremaju dinukleotidni indeksi
 * @param sequence Referenca na string u koji se sprema cijela bazna sekvenca
 * @param chr_number Broj kromosoma koji se učitava
 *
 * @throws runtime_error Ako se ulazna datoteka ne može otvoriti
 */
void load_chr_seq_to_dinuc_vector(vector<int>& O, string& s, int chr_number);


/**
 * @brief Ekstrahira kontinuirane CpG otoke iz vektora tvrdih stanja.
 *
 * Ulazni vektor stanja definira za svaki dinukleotid pripada li CpG
 * stanju (1) ili ne (0). Funkcija pronalazi maksimalne uzastopne segmente
 * CpG stanja i pretvara ih u CpG otoke izražene u baznim koordinatama.
 *
 * Mapiranje koordinata:
 *  - states[i] se odnosi na i-ti dinukleotid (0-based)
 *  - dinukleotid t (1-based) pokriva baze [t, t+1]
 *
 * @param islands Izlazni vektor CpG otoka (1-based bazne koordinate)
 * @param states Vektor tvrdih stanja po dinukleotidu
 *               (0 = non-CpG, 1 = CpG)
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
 * udaljeni najviše MERGE_DISTANCE baza, kako bi se izbacili prekratki otokci i
 * smanjio broj fragmentiranih otoka.
 * 
 * @param islands Referenca na vektor CpG otoka
 */
void filter_lenght_and_merge_close_islands(vector<CpgRegion>& islands);


/**
 * @brief Filtrira predviđene CpG otoke na temelju GC sadržaja i
 * omjera opaženih i očekivanih CpG dinukleotida (CpG O/E).
 *
 * Za svaki CpG otok računa se:
 *  - GC sadržaj = (broj C + broj G) / duljina otoka
 *  - CpG O/E = (broj opaženih "CG" dinukleotida * duljina) / (broj C * broj G)
 *
 * Otok se zadržava samo ako zadovoljava minimalne pragove
 * MIN_GC_CONTENT i MIN_CPG_OE.
 *
 * Koordinate CpG otoka su 1-based i uključive, dok je ulazna sekvenca
 * pohranjena kao 0-based string baza.
 *
 * @param sequence Bazna sekvenca kromosoma (A,C,G,T), 0-based indeksirana
 * @param islands Vektor CpG otoka u baznim koordinatama (modificira se in-place)
 */
void filter_by_content(const string& sequence, vector<CpgRegion>& islands);

