#pragma once


/**
 * Globalne konstante HMM-a
 *
 * Model ima:
 *  - 2 stanja:
 *      0 = Background
 *      1 = CpG island
 *
 *  - 16 emisijskih simbola koji predstavljaju dinukleotide:
 *      AA, AC, AG, AT,
 *      CA, CC, CG, CT,
 *      GA, GC, GG, GT,
 *      TA, TC, TG, TT
 */
constexpr int NSTATE = 2;
constexpr int NSYM   = 16;


/**
 * Struktura za pohranu HMM parametara
 *
 * A[i][j]  - prijelazna vjerojatnost iz stanja i u stanje j
 * B[i][k]  - emisijska vjerojatnost stanja i za za dinukleotid k (0..15)
 * pi[i]    - inicijalna vjerojatnost stanja i
 * chromosome označava kromosom na kojem je model zadnje treniran
 */
struct HMM {
    double A[NSTATE][NSTATE];
    double B[NSTATE][NSYM];
    double pi[NSTATE];
    int chromosome;
};


/**
 * Pretvara dinukleotid (prev, cur) u indeks 0..15:
 * AA=0, AC=1, AG=2, AT=3,
 * CA=4, CC=5, CG=6, CT=7,
 * GA=8, GC=9, GG=10, GT=11,
 * TA=12, TC=13, TG=14, TT=15
 *
 * Vraća -1 ako bilo koji znak nije A/C/G/T.
 */
inline int di_index(char a, char b) {
    int x =
        (a=='A')?0:(a=='C')?1:(a=='G')?2:(a=='T')?3:-1;
    int y =
        (b=='A')?0:(b=='C')?1:(b=='G')?2:(b=='T')?3:-1;
    return (x < 0 || y < 0) ? -1 : (x << 2) | y; // doslovno x*4+y
}


/**
 * CpgRegion struktura za pohranu koordinata CpG otoka
 * start      - početna pozicija (1-based)
 * end        - završna pozicija (1-based)
 * chromosome - kroj kromosoma
 */
struct CpgRegion { 
    int start;
    int end;
    int chromosome;
};


/**
 * Struktura za pohranu koordinata malih slova u genomu
 * start - početna pozicija (1-based)
 * end   - završna pozicija (1-based)
 * 
 * Napomena: stavljeno u 1-based radi konzistentnosti sa CpgRegion strukturama
 * koje su fiksno postavljene u 1-based sustav zbog USC koordinata.
 */
struct lowerCaseRegions {
    int start;
    int end;
};
