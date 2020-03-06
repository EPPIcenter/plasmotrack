//
// Created by Maxwell Murphy on 3/4/20.
//

#include "CombinationIndicesGenerator.h"

CombinationIndicesGenerator::CombinationIndicesGenerator(int n, int r) : completed(n < 1 || r > n), n_(n), r_(r) {
    for (int c = 0; c < r_; ++c)
        curr.push_back(c);
}

CombinationIndicesGenerator::combination_t CombinationIndicesGenerator::next() noexcept {
    combination_t ret = curr;
    completed = true;
    for (int i = r_ - 1; i >= 0; --i)
        if (curr[i] < n_ - r_ + i)
        {
            int j = curr[i] + 1;
            while (i <= r_)
                curr[i++] = j++;
            completed = false;
            generated++;
            break;
        }

    return ret;
}