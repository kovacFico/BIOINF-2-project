# CpG Island Detection using Hidden Markov Models

Ovaj projekt implementira cjelokupni **pipeline za detekciju CpG otoka** u ljudskom genomu korištenjem
**skrivenog Markovljevog modela (HMM)**.

Projekt obuhvaća:
- predobradu genomskih podataka
- inicijalizaciju HMM parametara
- treniranje modela Baum–Welch algoritmom
- dekodiranje sekvenci posteriornim pristupom
- evaluaciju predikcija u odnosu na poznate CpG otoke

Projekt je izrađen u sklopu kolegija **Bioinformatika 2**.

---

## 📁 Struktura projekta
```
.
├── src/
│ ├── main.cpp # Launcher aplikacije
│ │
│ ├── algorithms/
│ │ ├── baum_welch.cpp
│ │ ├── forward_backward.cpp
│ │
│ ├── apps/
│ │ ├── preprocess.cpp
│ │ ├── hmm_params_init.cpp
│ │ ├── train.cpp
│ │ └── decode_and_evaluation.cpp
│ │
│ ├── evaluation/
│ │ ├── evaluation.cpp
│ │
│ ├── hmm/
│ │ ├── hmm_io.cpp
│ │ ├── hmm.cpp
│ │
│ ├── postprocesing/
│ │ ├── decoded_postprocesing.cpp
│ │
│ ├── preprocesing/
│ │ ├── genome_preprocesing.cpp
│ │
│
├── include/
| ├── hmm/
| ├── evaluation/
| ├── algorithms/
| └── utils/
|
├── data/
│ └── ncbi_dataset/
|
├── output/
│ └── (generirane datoteke)
|
├── Makefile
└── README.md
```

---

## ⚙️ Preduvjeti

- C++17 kompatibilan kompajler (`g++`)
- Standardna C++ biblioteka
- Linux / macOS okruženje (projekt nije testiran na Windowsu)

---

## 🔧 Kompilacija

U korijenskom direktoriju projekta pokrenuti:
```
bash 
    make
```
Time se stvara izvršna datoteka:

bin/launcher


Za čišćenje build datoteka:

make clean

## ▶️ Pokretanje pipeline-a

Pokretanje cijelog pipeline-a (predobrada → treniranje → dekodiranje):

./bin/launcher


Launcher redom poziva sljedeće faze:

1. preprocess() – priprema genoma i CpG anotacija

2. hmm_params_init() – inicijalizacija HMM parametara

3. train_hmm() – treniranje HMM-a po kromosomima

4. decode_and_evaluate() – dekodiranje i evaluacija

## 🧠 Arhitektura pipeline-a

Pipeline je namjerno podijeljen u **više zasebnih izvršnih programa**
kako bi se izbjeglo prekomjerno korištenje memorije.

Svaka faza se izvršava u **posebnom procesu**, čime se:
- osigurava oslobađanje RAM-a nakon završetka faze
- izbjegava akumulacija memorije kod obrade velikih kromosoma
- omogućuje stabilno izvođenje na standardnim računalima

## 🧠 Metodologija

Model: 2-stanjni HMM (pozadina / CpG)

Emisije: dinukleotidi permutacijom A, C, G, T

Treniranje: Baum–Welch sa skaliranim forward/backward algoritmom

Dekodiranje: posteriorni pristup

Evaluacija: na razini CpG otoka i na razini parova baza

## 📌 Napomene

Putanje do ulaznih podataka trenutno su zadane u kodu stoga je
potrebna prilagodba ovisno gdje su postavljeni folderi. Po kodu
trenutno 'output' i 'data' folderi su izvan 'src' foldera u kojemu
je kod.

Projekt je optimiziran za velike sekvence (cijeli kromosomi)

Skaliranje se koristi radi numeričke stabilnosti

## 👤 Autor

Projekt izradio: Filip Barić, Filip Kovač

Studij: Računalstvo

Kolegij: Bioinformatika 2