#include "./forward_backward.hpp"    

double forward_scaled(
    const vector<int>& O,
    const HMM& hmm,
    vector<array<double, NSTATE>>& alpha,
    vector<double>& c
) {
    int T = (int)O.size();
    alpha.assign(T, {});
    c.assign(T, 0.0);

    // t = 0
    c[0] = 0.0;
    for (int i = 0; i < NSTATE; i++) {
        alpha[0][i] = hmm.pi[i] * hmm.B[i][O[0]];
        c[0] += alpha[0][i];
    }

    // AKO je suma 0 ili nije finite -> postavi uniformno i nastavi (sprječava 0*inf => NaN)
    if (!isfinite(c[0]) || c[0] <= 0.0) {
        for (int i = 0; i < NSTATE; i++) alpha[0][i] = 1.0 / NSTATE;
        c[0] = 1.0;
    } else {
        c[0] = 1.0 / c[0];
        if (!isfinite(c[0]) || c[0] > 1e300) c[0] = 1e300;
        for (int i = 0; i < NSTATE; i++) alpha[0][i] *= c[0];
    }

    // t = 1..T-1
    for (int t = 1; t < T; t++) {
        c[t] = 0.0;
        for (int j = 0; j < NSTATE; j++) {
            double sum = 0.0;
            for (int i = 0; i < NSTATE; i++)
                sum += alpha[t-1][i] * hmm.A[i][j];

            alpha[t][j] = sum * hmm.B[j][O[t]];
            c[t] += alpha[t][j];
        }

        // fallback ako suma pukne
        if (!isfinite(c[t]) || c[t] <= 0.0) {
            for (int j = 0; j < NSTATE; j++) alpha[t][j] = 1.0 / NSTATE;
            c[t] = 1.0;
        } else {
            c[t] = 1.0 / c[t];
            if (!isfinite(c[t]) || c[t] > 1e300) c[t] = 1e300;
            for (int j = 0; j < NSTATE; j++) alpha[t][j] *= c[t];
        }
    }

    double loglik = 0.0;
    for (int t = 0; t < T; t++) {
        // log(c[t]) je OK jer smo osigurali c[t] > 0 i finite
        loglik -= log(c[t]);
    }
    return loglik;
}

double forward_scaled_masked(
    const vector<int>& O,
    const HMM& hmm,
    const vector<array<double, NSTATE>>& state_mask,
    vector<array<double, NSTATE>>& alpha,
    vector<double>& c
) {
    int T = (int)O.size();
    alpha.assign(T, {});
    c.assign(T, 0.0);

    // t = 0
    c[0] = 0.0;
    for (int i = 0; i < NSTATE; i++) {
        alpha[0][i] = hmm.pi[i] * hmm.B[i][O[0]] * state_mask[0][i];
        c[0] += alpha[0][i];
    }

    if (!isfinite(c[0]) || c[0] <= 0.0) {
        double mask_sum = 0.0;
        for (int i = 0; i < NSTATE; i++) mask_sum += state_mask[0][i];
        if (mask_sum <= 0.0) {
            for (int i = 0; i < NSTATE; i++) alpha[0][i] = 1.0 / NSTATE;
        } else {
            for (int i = 0; i < NSTATE; i++) alpha[0][i] = state_mask[0][i] / mask_sum;
        }
        c[0] = 1.0;
    } else {
        c[0] = 1.0 / c[0];
        if (!isfinite(c[0]) || c[0] > 1e300) c[0] = 1e300;
        for (int i = 0; i < NSTATE; i++) alpha[0][i] *= c[0];
    }

    // t = 1..T-1
    for (int t = 1; t < T; t++) {
        c[t] = 0.0;
        for (int j = 0; j < NSTATE; j++) {
            double sum = 0.0;
            for (int i = 0; i < NSTATE; i++)
                sum += alpha[t-1][i] * hmm.A[i][j];

            alpha[t][j] = sum * hmm.B[j][O[t]] * state_mask[t][j];
            c[t] += alpha[t][j];
        }

        if (!isfinite(c[t]) || c[t] <= 0.0) {
            double mask_sum = 0.0;
            for (int j = 0; j < NSTATE; j++) mask_sum += state_mask[t][j];
            if (mask_sum <= 0.0) {
                for (int j = 0; j < NSTATE; j++) alpha[t][j] = 1.0 / NSTATE;
            } else {
                for (int j = 0; j < NSTATE; j++) alpha[t][j] = state_mask[t][j] / mask_sum;
            }
            c[t] = 1.0;
        } else {
            c[t] = 1.0 / c[t];
            if (!isfinite(c[t]) || c[t] > 1e300) c[t] = 1e300;
            for (int j = 0; j < NSTATE; j++) alpha[t][j] *= c[t];
        }
    }

    double loglik = 0.0;
    for (int t = 0; t < T; t++) {
        loglik -= log(c[t]);
    }
    return loglik;
}


void backward_scaled(
    const vector<int>& O,
    const HMM& hmm,
    const vector<double>& c,
    vector<array<double, NSTATE>>& beta
) {
    int T = O.size();
    beta.assign(T, {});

    for (int i = 0; i < NSTATE; i++)
        beta[T-1][i] = c[T-1];

    for (int t = T-2; t >= 0; t--) {
        for (int i = 0; i < NSTATE; i++) {
            beta[t][i] = 0.0;
            for (int j = 0; j < NSTATE; j++)
                beta[t][i] += hmm.A[i][j] * hmm.B[j][O[t+1]] * beta[t+1][j];
            beta[t][i] *= c[t];
        }
    }
}


void backward_scaled_masked(
    const vector<int>& O,
    const HMM& hmm,
    const vector<array<double, NSTATE>>& state_mask,
    const vector<double>& c,
    vector<array<double, NSTATE>>& beta
) {
    int T = O.size();
    beta.assign(T, {});

    for (int i = 0; i < NSTATE; i++)
        beta[T-1][i] = c[T-1] * state_mask[T-1][i];

    for (int t = T-2; t >= 0; t--) {
        for (int i = 0; i < NSTATE; i++) {
            beta[t][i] = 0.0;
            for (int j = 0; j < NSTATE; j++)
                beta[t][i] += hmm.A[i][j] * hmm.B[j][O[t+1]] * state_mask[t+1][j] * beta[t+1][j];
            beta[t][i] *= c[t];
        }
    }
}
