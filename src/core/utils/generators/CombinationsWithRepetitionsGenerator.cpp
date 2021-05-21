//
// Created by Maxwell Murphy on 4/29/20.
//

#include "CombinationsWithRepetitionsGenerator.h"


namespace transmission_nets::core::utils::generators {
    CombinationsWithRepetitionsGenerator::CombinationsWithRepetitionsGenerator(int nChoices, int k) : nChoices_(nChoices),
                                                                                                      k_(k) {
        curr = std::vector<int>(k_, 0);
        curr[0] = -1;
    }


    void CombinationsWithRepetitionsGenerator::reset(int nChoices, int k) noexcept {
        k_ = k;
        nChoices_ = nChoices;
        std::fill(curr.begin(), curr.end(), 0);
        curr.resize(k_, 0);
        curr[0] = -1;
        generated = 0;
    }


    void CombinationsWithRepetitionsGenerator::next() noexcept {
        generated++;

        int currIdx = 0;

        // fast break in the common case
        if (curr[currIdx] < nChoices_ - 1) {
            curr[currIdx]++;
            if (nChoices_  == 1) {
                completed = true;
            }

            if (curr[currIdx] == nChoices_ - 1 and k_ == 1) {
                completed = true;
            }

            return;
        }

        currIdx++;
        bool carry = true;
        while (carry) {
            if (curr[currIdx] == nChoices_ - 1) {
                currIdx++;
                continue;
            }

            carry = false;
            curr[currIdx]++;
            if(currIdx == k_ - 1 and curr[currIdx] == nChoices_ - 1) {
                completed = true;
                break;
            }

            int upper = curr[currIdx];
            currIdx--;
            while(currIdx >= 0) {
                curr[currIdx] = upper;
                currIdx--;
            }
        }
    }

    CombinationsWithRepetitionsGenerator::CombinationsWithRepetitionsGenerator() {
        nChoices_ = 1;
        k_ = 1;
        curr = std::vector<int>(k_, 0);
        curr[0] = -1;
    }


}



