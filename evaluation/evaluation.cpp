#include "./evaluation.hpp"


void island_based_evaluation(const vector<CpgRegion>& predicted_in, const vector<CpgRegion>& truth_in) {
    vector<CpgRegion> predicted = predicted_in;
    vector<CpgRegion> truth     = truth_in;

    sort(predicted.begin(), predicted.end(),
         [](const CpgRegion& a, const CpgRegion& b){ return a.start < b.start; });

    sort(truth.begin(), truth.end(),
         [](const CpgRegion& a, const CpgRegion& b){ return a.start < b.start; });

    int TP = 0;
    size_t i = 0; // predicted
    size_t j = 0; // truth

    while (i < predicted.size() && j < truth.size()) {
        const auto& p = predicted[i];
        const auto& t = truth[j];

        if (p.end < t.start) {
            i++;
            continue;
        }

        if (t.end < p.start) {
            j++;
            continue;
        }

        TP++;
        i++;
        j++;
    }

    int FP = (int)predicted.size() - TP;
    int FN = (int)truth.size()     - TP;

    double precision = 0.0;
    double recall    = 0.0;

    if (TP + FP > 0) precision = TP / double(TP + FP);
    if (TP + FN > 0) recall    = TP / double(TP + FN);

    cout << "=== Island-based evaluation ===\n";
    cout << "True Positive = " << TP << "\n";
    cout << "False Positive = " << FP << "\n";
    cout << "False Negative = " << FN << "\n";
    cout << "Precision = " << precision << "\n";
    cout << "Recall = " << recall << "\n";
}


void base_pair_evaluation(const vector<CpgRegion>& predicted, const vector<CpgRegion>& truth) {
    long long pred_len = 0;
    long long truth_len = 0;

    for (const auto& p : predicted)
        pred_len += (p.end - p.start + 1);

    for (const auto& t : truth)
        truth_len += (t.end - t.start + 1);

    long long overlap_len = 0;

    size_t i = 0, j = 0;
    while (i < predicted.size() && j < truth.size()) {
        int s = max(predicted[i].start, truth[j].start);
        int e = min(predicted[i].end,   truth[j].end);

        if (s <= e)
            overlap_len += (e - s + 1);

        if (predicted[i].end < truth[j].end)
            i++;
        else
            j++;
    }

    long long TP = overlap_len;
    long long FP = pred_len  - overlap_len;
    long long FN = truth_len - overlap_len;

    double precision = TP / double(TP + FP);
    double recall    = TP / double(TP + FN);

    cout << "=== Base-pair evaluation ===\n";
    cout << "True Positive = " << TP << "\n";
    cout << "False Positive = " << FP << "\n";
    cout << "False Negative = " << FN << "\n";
    cout << "Precision = " << precision << "\n";
    cout << "Recall = " << recall << "\n";
}