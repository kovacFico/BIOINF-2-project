#include <cmath>

#include "./baum_welch.hpp"

double baum_welch_iteration_multi_masked(
    const vector<vector<int>>& sequences,
    const vector<vector<array<double, NSTATE>>>& state_masks,
    HMM& hmm,
    double& ll
) {
    double A_num[NSTATE][NSTATE] = {{0}};
    double A_den[NSTATE] = {0};

    double B_num[NSTATE][NSYM] = {{0}};
    double B_den[NSTATE] = {0};

    int used_sequences = 0;

    for (size_t sidx = 0; sidx < sequences.size(); sidx++) {
        const auto& O = sequences[sidx];
        const auto& mask = state_masks[sidx];
        int T = (int)O.size();
        if (T < 2) continue;

        vector<array<double, NSTATE>> alpha, beta;
        vector<double> c;

        double lseq = forward_scaled_masked(O, hmm, mask, alpha, c);
        if (!isfinite(lseq)) continue;
        ll += lseq;

        backward_scaled_masked(O, hmm, mask, c, beta);

        for (int t = 0; t < T - 1; t++) {
            double xi_den = 0.0;
            for (int i = 0; i < NSTATE; i++)
                for (int j = 0; j < NSTATE; j++)
                    xi_den += alpha[t][i] * hmm.A[i][j] * hmm.B[j][O[t+1]] * mask[t+1][j] * beta[t+1][j];

            double gamma_den = 0.0;
            for (int i = 0; i < NSTATE; i++)
                gamma_den += alpha[t][i] * beta[t][i];

            if (!isfinite(xi_den) || xi_den <= 0.0) xi_den = 1e-300;
            if (!isfinite(gamma_den) || gamma_den <= 0.0) gamma_den = 1e-300;

            for (int i = 0; i < NSTATE; i++) {
                double gamma = (alpha[t][i] * beta[t][i]) / gamma_den;
                if (!isfinite(gamma) || gamma < 0.0) continue;

                A_den[i] += gamma;
                B_den[i] += gamma;
                B_num[i][O[t]] += gamma;

                for (int j = 0; j < NSTATE; j++) {
                    double xi = alpha[t][i] * hmm.A[i][j] * hmm.B[j][O[t+1]] * mask[t+1][j] * beta[t+1][j] / xi_den;
                    if (!isfinite(xi) || xi < 0.0) continue;
                    A_num[i][j] += xi;
                }
            }
        }

        {
            int t = T - 1;
            double gamma_den = 0.0;
            for (int i = 0; i < NSTATE; i++) gamma_den += alpha[t][i] * beta[t][i];
            if (!isfinite(gamma_den) || gamma_den <= 0.0) gamma_den = 1e-300;

            for (int i = 0; i < NSTATE; i++) {
                double gamma = (alpha[t][i] * beta[t][i]) / gamma_den;
                if (!isfinite(gamma) || gamma < 0.0) continue;

                B_den[i] += gamma;
                B_num[i][O[t]] += gamma;
            }
        }

        used_sequences++;
    }

    if (used_sequences == 0) {
        return ll;
    }

    const double A_PSEUDO = 1e-3;
    const double B_PSEUDO = 1e-2;
    const double B_FLOOR = 1e-6;

    for (int i = 0; i < NSTATE; i++) {
        if (!isfinite(A_den[i]) || A_den[i] <= 0.0) A_den[i] = 1e-300;
        if (!isfinite(B_den[i]) || B_den[i] <= 0.0) B_den[i] = 1e-300;

        double Aden = A_den[i] + A_PSEUDO * NSTATE;
        double Bden = B_den[i] + B_PSEUDO * NSYM;

        for (int j = 0; j < NSTATE; j++) {
            double num = A_num[i][j];
            if (!isfinite(num) || num < 0.0) num = 0.0;
            hmm.A[i][j] = (num + A_PSEUDO) / Aden;
        }

        double bsum = 0.0;
        for (int k = 0; k < NSYM; k++) {
            double num = B_num[i][k];
            if (!isfinite(num) || num < 0.0) num = 0.0;
            hmm.B[i][k] = (num + B_PSEUDO) / Bden;
            if (hmm.B[i][k] < B_FLOOR) hmm.B[i][k] = B_FLOOR;
            bsum += hmm.B[i][k];
        }
        if (bsum > 0.0) {
            for (int k = 0; k < NSYM; k++) {
                hmm.B[i][k] /= bsum;
            }
        }
    }

    return ll;
}

