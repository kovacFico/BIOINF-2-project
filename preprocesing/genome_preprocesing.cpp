#include "./genome_preprocesing.hpp"

vector<string> load_positive_cpg(const string &filename, vector<CpgRegion> &coords) {
    ifstream file(filename);
    if (!file) {
        cerr << "Ne mogu otvoriti pozitivan dataset!" << endl;
        exit(1);
    }

    vector<string> cpg_list;
    string line;
    string current_seq;
    int chromosome, start, end;

    while (getline(file, line)) {
        if (line.size() == 0) continue;

        if (line[0] == '>') {
            // ako vec imamo sekvencu, spremamo ju
            if (!current_seq.empty()) {
                cpg_list.push_back(current_seq);
                current_seq.clear();
            }

            size_t pos = line.find("range=chr");
            if (pos != string::npos) {
                pos += string("range=chr").size();
                size_t colon = line.find(':', pos);
                try {
                    chromosome = stoi(line.substr(pos, colon - pos));
                } catch (const invalid_argument &e) {
                    // naišli smo na zapis kromosoma koji nije 
                    // broj; prekidamo petlju jer su ostali samo
                    // X i Y
                    continue;
                }

                pos = colon + 1;
                size_t dash = line.find('-', pos);
                size_t space = line.find(' ', dash);
                start = stoll(line.substr(pos, dash - pos));
                end   = stoll(line.substr(dash + 1, space - dash - 1));

                coords.push_back({start, end, chromosome});
            }

        } else {
            current_seq += line;
        }
    }

    if (!current_seq.empty())
        cpg_list.push_back(current_seq);

    return cpg_list;
}


string load_chromosome(const string &filename, vector<lowerCaseRegions> &lowercaseCoords, int chr_number) {
    ifstream file(filename);
    if (!file) {
        cerr << "Ne mogu otvoriti kromosom!" << endl;
        exit(1);
    }

    string line, genome;
    bool in_correct_chr = false;
    bool in_lowercase = false;
    int start = -1;
    int pos = 1;

    while (getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (in_correct_chr) break;
            if (line.find("chromosome " + to_string(chr_number)) != string::npos) in_correct_chr = true;
        } else if (in_correct_chr) {
            for (char c : line) {
                if (isupper(c)) {
                    genome.push_back(c);

                    if (in_lowercase) {
                        lowercaseCoords.push_back({start, pos - 1});
                        in_lowercase = false;
                    }

                } else {
                    if (!in_lowercase) {
                        start = pos;
                        in_lowercase = true;
                    }
                }
                pos++;
            }            
        }
    }

    if (in_lowercase) lowercaseCoords.push_back({start, pos - 1});
    
    return genome;
}


string load_background(const string &filename, const vector<CpgRegion> &coords) {
    ifstream file(filename);
    if (!file) {
        cerr << "Ne mogu otvoriti kromosom!" << endl;
        exit(1);
    }

    string line, genome;
    while (getline(file, line)) {
        if (line.size() == 0) continue;
        if (line[0] == '>') {
            size_t pos = line.find("chromosome ");
            if (pos != string::npos) {
                pos += string("chromosome ").size();
                try {
                    stoi(line.substr(pos, pos + 1));
                    continue;
                } catch (const invalid_argument &e) {
                    // naišli smo na zapis kromosoma koji nije 
                    // broj; prekidamo petlju jer su ostali samo
                    // X i Y
                    break;
                }
            }
        }

        for (char c : line) {
            if (isupper(c)) genome.push_back(c);
        }
    }

    string cleaned;
    cleaned.reserve(genome.size());
    
    for (const auto &r : coords) {
        //PAZI: UCSC koordinate su 1-based, C++ je 0-based
        for (int i = r.start - 1; i <= r.end - 1 && i < (int)genome.size(); i++) {
            genome[i] = '\0';  // Postavljamo CpG otok na '\0' kako bi ga preskočili
        }
    }

    for (char c : genome) {
        if (c != '\0') cleaned.push_back(c);
    }
    return cleaned;
}


void open_output_files(
    int num_chromosomes, 
    const string &output_dir, 
    vector<ofstream> &chromosome_out_files, 
    ofstream &out1, 
    ofstream &out2, 
    ofstream &coords_out
) {

    chromosome_out_files.resize(num_chromosomes);
    
    for (int chr = 1; chr <= num_chromosomes; chr++) {
        string chr_filename = output_dir + "/" + to_string(chr) + (chr <= 16 ? "_train_chr.txt" : "_test_chr.txt");
        chromosome_out_files[chr - 1].open(chr_filename);
        if (!chromosome_out_files[chr - 1]) {
            cerr << "Greška pri otvaranju fajla za kromosom " << chr << endl;
            exit(1);
        }
    }

    out1.open(output_dir + "/clean_positive.txt");
    out2.open(output_dir + "/clean_background.txt");
    coords_out.open(output_dir + "/coords.txt");

    if (!out1 || !out2 || !coords_out) {
        cerr << "Greška pri otvaranju izlaznih fajlova!" << endl;
        exit(1);
    }
}