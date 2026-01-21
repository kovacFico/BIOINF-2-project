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


/**
 * @brief Ucitava sekvencu + lowercase intervale iz *_train_chr fajla.
 * 
 * @param filename Putanja do fajla.
 * @param seq Referenca na string u koji se učitava sekvenca.
 * @param lc Referenca na vektor lowercase regija koje se učitavaju.
 */
void get_chromosome_and_lowercase_regions(const string& filename, string& seq, vector<lowerCaseRegions>& lc);


/**
 * @brief Pretvara sekvencu kromosoma u dinukleotidna opažanja i sprema
 * ih kao vrijednosti [1, 16] u vektor O. Pogledaj funkciju di_index u 
 * structs_consts_functions.hpp za mapiranje dinukleotida na indekse.
 * 
 * @param s Sekvenca kromosoma (1-based indeksi).
 * @param start Početna pozicija (1-based).
 * @param end Završna pozicija (1-based).
 * 
 * @return Vektor dinukleotidnih opažanja.
 */
vector<int> seq_to_dinuc(const string& s, int start, int end);


/**
 * @brief Mapira original 1-based start/end coord -> komprimirana 1-based start/end coord (uppercase-only)
 * 
 * @param lc Vektor lowercase regija.
 * @param orig_pos Originalna 1-based pozicija.
 * 
 * @return Komprimirana 1-based pozicija.
 */
int map_orig_to_comp(const vector<lowerCaseRegions>& lc, int orig_pos);


/**
 * Mapira originalne koordinate CpG otoka na koordinate u komprimiranoj sekvenci
 * (bez lowercase regija). Odbija od orginalnih koordinata dijelove koji su u 
 * lowercase regijama, kako bi se dobile koordinate u komprimiranoj sekvenci.
 *
 * @param seq Komprimirana sekvenca kromosoma (samo uppercase).
 * @param lc Vektor lowercase regija.
 * @param orig_coords Vektor originalnih koordinata CpG otoka.
 * 
 * @return Vektor koordinata CpG otoka u komprimiranoj sekvenci.
 */
vector<CpgRegion> map_orig_coords_to_compressed(
    const string& seq,
    const vector<lowerCaseRegions>& lc,
    const vector<CpgRegion>& orig_coords
);


/**
 * Izgrađuje trening sekvence dinukleotida i pripadajuće maske dozvoljenih HMM stanja
 * za semi-supervizirano treniranje CpG HMM-a.
 *
 * Funkcija kombinira:
 *  - DNA sekvencu kromosoma (uppercase-only),
 *  - komprimirane koordinate CpG regija,
 *  - lokalni kontekst oko CpG regija (±NEG_MARGIN),
 * kako bi za svaki dinukleotid odredila koja su HMM stanja dozvoljena.
 *
 * Ideja semi-supervizije:
 *  - dinukleotidi koji se nalaze unutar poznatih CpG regija su "clamped" na CpG stanje,
 *  - dinukleotidi koji su dovoljno daleko od svih CpG regija su "clamped" na non-CpG stanje,
 *  - dinukleotidi u prijelaznim zonama (blizu granica CpG regija) ostaju neoznačeni,
 *    tj. dopuštaju oba stanja kako bi HMM sam naučio prijelaze.
 *
 * Zbog memorijskih i računalnih ograničenja, sekvenca se dijeli u nezavisne
 * chunkove fiksne maksimalne duljine (izražene u broju dinukleotida).
 * Svaki chunk se tretira kao zasebna trening sekvenca u Baum–Welch algoritmu.
 *
 * @param s Uppercase DNA sekvenca kromosoma (1-based indeksiranje se koristi logički).
 * @param coords_chr Koordinate poznatih CpG regija, mapirane u komprimirani prostor.
 * @param sequences Izlazni vektor dinukleotidnih opažanja (jedan vektor po chunku).
 * @param masks Izlazni vektor maski dozvoljenih stanja (paralelan s `sequences`).
 */
void build_masked_sequences(
    const string& s,
    const vector<CpgRegion>& coords_chr,
    vector<vector<int>>& sequences,
    vector<vector<array<double, NSTATE>>>& masks
);