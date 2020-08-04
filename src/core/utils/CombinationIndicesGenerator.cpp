//
// Created by Maxwell Murphy on 3/4/20.
//

#include <numeric>
#include "CombinationIndicesGenerator.h"

CombinationIndicesGenerator::CombinationIndicesGenerator(int n, int r) : completed(n < 1 or r > n or r == 0), n_(n),
                                                                         r_(r) {
    curr.resize(r_);
    std::iota(curr.begin(), curr.end(), 0);
//    curr.clear();
//    for (int i = 0; i < r_; ++i) {
//        curr.push_back(i);
//    }
}

void CombinationIndicesGenerator::reset(int n, int r) {
    completed = n < 1 or r > n or r == 0;
    generated = 1;

    n_ = n;
    r_ = r;

    curr.resize(r_);
    std::iota(curr.begin(), curr.end(), 0);
//    curr.clear();
//    for (int i = 0; i < r_; ++i) {
//        curr.push_back(i);
//    }
}

void CombinationIndicesGenerator::next() noexcept {
    assert(!completed);

    completed = true;
    for (int i = r_ - 1; i >= 0; --i)
        if (curr.at(i) < n_ - r_ + i) {
            int j = curr.at(i) + 1;
            while (i < r_) {
                curr.at(i++) = j++;
            }
            completed = false;
            generated++;
            break;
        }
}

CombinationIndicesGenerator::CombinationIndicesGenerator() {
    completed = true;
    n_ = 0;
    r_ = 0;
}

//CombinationIndicesGenerator::combination_t CombinationIndicesGenerator::curr = {};