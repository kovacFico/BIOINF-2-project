#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "../utils/structs_consts_functions.hpp"

using namespace std;


/**
 * Učitavanje pozitivnih CpG otoka iz FASTA fajla. Također parsira koordinate iz header linija
 * Rezultat: vektor sekvenci CpG otoka i popunjen vektor koordinata coords.
 * 
 * @param filename - ime FASTA fajla
 * @param coords - referenca na vektor za pohranu koordinata
 * 
 * @return vektor sekvenci CpG otoka
 */
vector<string> load_positive_cpg(const string &filename, vector<CpgRegion> &coords);


/**
 * Učitava kromosom iz FASTA fajla, uključujući samo velika slova i ne čisti CpG otoke.
 * Također sprema raspone malih slova u genomu. Kromosome kasnije koristimo u Baum-Welch treningu.
 * 
 * @param filename Ime FASTA fajla
 * @param lowercaseCoords Referenca na vektor za pohranu koordinata malih slova
 * @param chr_number Broj kromosoma za učitavanje
 * 
 * @return string Sekvenca kromosoma
 */
string load_chromosome(const string &filename, vector<lowerCaseRegions> &lowercaseCoords, int chr_number);


/**
 * Učitava genom oćišćen od CpG otoka i malih slova iz .txt doteteke.
 * Baziran je na koordinatama CpG otoka. Koristimo ga za inicijalizaciju početnog stanja HMM-a.
 * 
 * @param filename Ime doteteke
 * @param coords Vektor koordinata CpG otoka
 * 
 * @return string Sekvenca pozadinskog genoma
 * 
 * Napomena: Pretpostavlja se da je sekvenca već očišćena i spremna za upotrebu.
 */
string load_background(const string &filename, const vector<CpgRegion> &coords);


/**
 * Otvara izlazne fajlove za svaki kromosom, pozitivne CpG otoke, pozadinski genom i koordinate
 * 
 * @param num_chromosomes Broj kromosoma
 * @param output_dir Direktorij za izlazne fajlove
 * @param chromosome_out_files Referenca na vektor ofstream objekata za kromosome
 * @param out1 Referenca na ofstream za pozitivne CpG otoke
 * @param out2 Referenca na ofstream za pozadinski genom
 * @param coords_out Referenca na ofstream za koordinate CpG otoka
 */
void open_output_files(
    int num_chromosomes, 
    const string &output_dir, 
    vector<ofstream> &chromosome_out_files, 
    ofstream &out1, 
    ofstream &out2, 
    ofstream &coords_out
);