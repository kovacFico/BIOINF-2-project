#include "./decoded_postprocesing.hpp"


void keep_and_clip(vector<CpgRegion>& islands, int keep_start_bp, int keep_end_bp) {
    vector<CpgRegion> out;
    out.reserve(islands.size());

    for (const auto r : islands) {
        int s = max(r.start, keep_start_bp);
        int e = min(r.end, keep_end_bp);
        if (e >= s) out.push_back({s, e, r.chromosome});
    }

    islands.swap(out);
}


void trim_islands_with_posterior(
    vector<CpgRegion>& islands,
    const vector<double>& posterior_c,
    int base_shift,
    double trim_th
) {
    vector<CpgRegion> trimmed;
    trimmed.reserve(islands.size());

    const int max_d = static_cast<int>(posterior_c.size());

    for (const auto& r : islands) {
        int d_start = max(1, r.start - base_shift);
        int d_end   = min(max_d, r.end - base_shift - 1);

        // Trim strana dok je posterior ispod praga
        while (d_start <= d_end && posterior_c[d_start - 1] < trim_th) d_start++;
        while (d_end >= d_start && posterior_c[d_end - 1] < trim_th) d_end--;

        if (d_start <= d_end) {
            trimmed.push_back({
                d_start + base_shift, 
                d_end + 1 + base_shift, 
                r.chromosome
            });
        }
    }

    islands.swap(trimmed);
}


void load_chr_seq_to_dinuc_vector(vector<int>& O, string& s, int chr_number) {
    ifstream in("../output/" + to_string(chr_number) + "_test_chr.txt");
    if (!in) {
        cerr << "Ne mogu otvoriti test fajl za kromosom " << chr_number << "\n";
        exit(1);
    }
    getline(in, s);

    O.reserve(s.size() - 1);
    for (size_t i = 1; i < s.size(); i++) {
        int x = di_index(s[i - 1], s[i]);
        if (x != -1) O.push_back(x);
    }
}


void extract_cpg_islands(vector<CpgRegion>& islands, vector<int>& states) {
    int start_d = -1;

    for (int i = 0; i < (int)states.size(); i++) {
        if (states[i] == 1 && start_d == -1) start_d = i + 1;

        if ((states[i] == 0 || i == (int)states.size() - 1) && start_d != -1) {
            int end_d = (states[i] == 0) ? i : i + 1; // dinukleotid end (1-based)
            int start_bp = start_d;
            int end_bp = end_d + 1;

            islands.push_back({start_bp, end_bp, 0});
            start_d = -1;
        }
    }
}


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



