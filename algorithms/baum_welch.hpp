#pragma once

#include "../utils/structs_consts_functions.hpp"
#include "../algorithms/forward_backward.hpp"

using namespace std;

double baum_welch_iteration_multi_masked(
    const vector<vector<int>>& sequences,
    const vector<vector<array<double, NSTATE>>>& state_masks,
    HMM& hmm,
    double& ll
);
