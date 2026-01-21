#include "./decoded_postprocesing.hpp"

// Pretvori binarnu sekvencu stanja (0/1 po dinukleotidu) u CpG intervale u baznim koordinatama.
void extract_cpg_islands(vector<CpgRegion>& islands, vector<int>& states) {
    int start_d = -1; // start u dinukleotid indeksima (1-based)

    for (int i = 0; i < (int)states.size(); i++) {
        if (states[i] == 1 && start_d == -1) start_d = i + 1;

        if ((states[i] == 0 || i == (int)states.size() - 1) && start_d != -1) {
            int end_d = (states[i] == 0) ? i : i + 1; // dinukleotid end (1-based)

            // MAPIRANJE: dinukleotid [start_d, end_d] pokriva baze [start_d, end_d+1]
            int start_bp = start_d;
            int end_bp = end_d + 1;

            islands.push_back({start_bp, end_bp, 0}); // chr nije bitan ovdje
            start_d = -1;
        }
    }
}

// Vrati predikcije iz uppercase-only sekvence natrag u originalne koordinate dodavanjem lowercase offseta.
void move_predicted_based_on_lowercase(vector<CpgRegion>& predicted, int chr_number) {
    vector<CpgRegion> lowercaseCoords;
    ifstream in("../output/" + to_string(chr_number) + "_test_chr.txt");
    string s;
    getline(in, s);

    while(getline(in, s)) {
        size_t space = s.find(' ', 0);
        int start = stoll(s.substr(0, space));
        int end   = stoll(s.substr(space + 1, s.size() - space - 1));
        lowercaseCoords.push_back({start, end});
    }


    for (auto& p : predicted) {
        int offset = 0;

        for (const auto& lc : lowercaseCoords) {
            if (lc.start <= p.start + offset) {
                offset += (lc.end - lc.start + 1);
            } else {
                break;
            }
        }

        p.start += offset;
        p.end += offset;

    }
}

// Spoji bliske CpG otoke i odbaci prekratke.
void filter_lenght_and_merge_close_islands(vector<CpgRegion>& islands) {
    if (islands.empty()) return;

    sort(islands.begin(), islands.end(), 
        [](const CpgRegion& a, const CpgRegion& b) {
            return a.start < b.start;
        });

    vector<CpgRegion> temp;
    temp.push_back(islands[0]);

    for (int i = 1; i < islands.size(); i++) {
        if (islands[i].start - temp.back().end <= MERGE_DISTANCE) {
            temp.back().end = max(temp.back().end, islands[i].end);
        } else {
            temp.push_back(islands[i]);
        }
    }

    islands.clear();

    for (const auto& region : temp) {
        if (region.end - region.start + 1 >= MIN_CPG_LEN) {
            islands.push_back(region);
        }
    }
}

// Biološki filter: zadrži islands koji zadovoljavaju GC content i CpG O/E pragove.
void filter_by_content(const string& sequence, vector<CpgRegion>& islands) {
    if (islands.empty()) return;

    vector<CpgRegion> filtered;
    filtered.reserve(islands.size());

    for (const auto& region : islands) {
        int start = max(1, region.start);
        int end = min((int)sequence.size(), region.end);
        int len = end - start + 1;
        if (len <= 0) continue;

        int count_c = 0;
        int count_g = 0;
        int count_cg = 0;

        char prev = '\0';
        for (int pos = start; pos <= end; pos++) {
            char base = toupper(sequence[(size_t)(pos - 1)]);
            if (base == 'C') count_c++;
            if (base == 'G') count_g++;
            if (prev == 'C' && base == 'G') count_cg++;
            prev = base;
        }

        double gc_content = (count_c + count_g) / double(len);
        double oe = 0.0;
        if (count_c > 0 && count_g > 0) {
            oe = (count_cg * double(len)) / (count_c * double(count_g));
        }

        if (gc_content >= MIN_GC_CONTENT && oe >= MIN_CPG_OE) {
            filtered.push_back(region);
        }
    }

    islands.swap(filtered);
}

void load_true_islands(vector<CpgRegion>& islands, int chr_number) {
    ifstream in("../output/coords.txt");
    int chr, start, end;

    while (in >> chr >> start >> end) {
        if (chr == chr_number) {
            islands.push_back({start, end});
        }
    } 
}



