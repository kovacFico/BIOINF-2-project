#pragma once

#include <string>

#include "./forward_backward.hpp"

using namespace std;


/**
 * @brief Izračunava posteriorne vjerojatnosti P(Z_t = CpG | Oseg) za svaki dinukleotid da je CpG stanje 
 * koristeći skalirani forward i backward algoritam. 
 * 
 * @param Oseg Vektor opažanja (dinukleotidi) za segment
 * @param hmm HMM model s treniranim parametrima
 * 
 * @return vector<double> Vektor posteriornih vjerojatnosti za svaki dinukleotid u segmentu
 * 
 * Napomena: Funkcija koristi skalirane verzije forward i backward algoritama kako bi se izbjegle numeričke nestabilnosti.
 */
vector<double> compute_posterior_c(const vector<int>& Oseg, const HMM& hmm);


/**
 * @brief Pretvara posteriorne vjerojatnosti u tvrda stanja (0/1) koristeći histerezis
 * s definiranim pragovima ulaska i izlaska. Ovo pomaže u smanjenju "treperenja"
 * predviđenih CpG otoka.
 * 
 * @param posterior_c Vektor posteriornih vjerojatnosti za svaki dinukleotid
 * @param enter_th Prag ulaska u CpG stanje
 * @param exit_th Prag izlaska iz CpG stanje
 * 
 * @return vector<int> Vektor tvrdih stanja (0 = non-CpG, 1 = CpG) za svaki dinukleotid
 */
vector<int> decode_hysteresis(const vector<double>& posterior_c, double enter_th, double exit_th);


/**
 * @brief Obrada jednog preklapajućeg prozora dinukleotida radi predikcije CpG otoka.
 *
 * Funkcija:
 *  - izdvaja segment opažanja
 *  - računa posteriorne vjerojatnosti (forward/backward)
 *  - dekodira tvrda stanja s histerezom
 *  - ekstrahira CpG otoke
 *  - mapira ih u globalne bazne koordinate
 *  - uklanja rubne artefakte, trimma po posterioru i filtrira po sadržaju
 *
 * @param O Globalni vektor dinukleotidnih opažanja
 * @param sequence Globalna bazna sekvenca kromosoma
 * @param hmm Trenirani HMM model
 * @param start_d Početni indeks dinukleotida prozora (0-based)
 * @param end_d Završni indeks dinukleotida prozora (exclusive)
 * @param T Ukupan broj dinukleotida u sekvenci
 * @param POST_ENTER Prag ulaska u CpG stanje (histerezis)
 * @param POST_EXIT Prag izlaska iz CpG stanja (histerezis)
 * @param POST_TRIM Prag posteriora za trimanje rubova CpG otoka
 * @param OVERLAP Broj dinukleotida preklapanja između susjednih prozora
 *
 * @return vector<CpgRegion> Lista predviđenih CpG otoka u globalnim baznim koordinatama
 */
vector<CpgRegion> process_window(
    const vector<int>& O,
    const string& sequence,
    const HMM& hmm,
    int start_d,
    int end_d,
    int T,
    double POST_ENTER, 
    double POST_EXIT, 
    double POST_TRIM, 
    int OVERLAP
);
