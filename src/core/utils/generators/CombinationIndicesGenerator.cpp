//
// Created by Maxwell Murphy on 3/4/20.
//

#include "CombinationIndicesGenerator.h"

#include <cassert>
#include <numeric>

namespace transmission_nets::core::utils::generators {
    CombinationIndicesGenerator::CombinationIndicesGenerator(std::size_t n, std::size_t r) : completed(n < 1 or r > n or r == 0),
                                                                                             n_(n),
                                                                                             r_(r) {
        curr.resize(r_);
        std::iota(curr.begin(), curr.end(), 0);
    }

    void CombinationIndicesGenerator::reset(std::size_t n, std::size_t r) {
        completed = n < 1 or r > n or r == 0;
        generated = 1;

        n_ = n;
        r_ = r;

        curr.resize(r_);
        std::iota(curr.begin(), curr.end(), 0);
    }

    void CombinationIndicesGenerator::reset() {
        reset(n_, r_);
    }

    void CombinationIndicesGenerator::next() noexcept {
        assert(!completed);
        completed = true;
        for (long i = (signed) r_ - 1; i >= 0; --i) {
            if (curr.at(i) < n_ - r_ + i) {
                unsigned int j = curr.at(i) + 1;
                while (i < (signed) r_) {
                    curr.at(i++) = j++;
                }
                completed = false;
                generated++;
                break;
            }
        }
    }

    void CombinationIndicesGenerator::advance(int n) noexcept {
        for (int i = 0; i < n; ++i) {
            next();
        }
    }

    CombinationIndicesGenerator::CombinationIndicesGenerator() {
        completed = true;
        n_        = 0;
        r_        = 0;
    }
}// namespace transmission_nets::core::utils::generators
