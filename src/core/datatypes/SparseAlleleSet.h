//
// Created by mmurphy on 5/28/24.
//

#ifndef SPARSEALLELESET_H
#define SPARSEALLELESET_H

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <cstring>

namespace transmission_nets::core::datatypes {
    class SparseAlleleSet {
        static constexpr unsigned int MAX_SIZE = 30;
    public:
        explicit SparseAlleleSet(const std::string& bitstr);
        explicit SparseAlleleSet(unsigned int totalAlleles);
        constexpr SparseAlleleSet(unsigned int totalAlleles, const std::array<unsigned char, MAX_SIZE>& set_positions) : total_alleles_(totalAlleles) {
            set_count_ = set_positions.size();
            std::memcpy(set_positions_, set_positions.data(), sizeof(set_positions_));
        }
        SparseAlleleSet(const SparseAlleleSet& other);
        SparseAlleleSet(SparseAlleleSet&& other) noexcept;
        SparseAlleleSet();

        SparseAlleleSet& operator=(const SparseAlleleSet& other);
        SparseAlleleSet& operator=(SparseAlleleSet&& other) noexcept ;

        friend std::ostream& operator<<(std::ostream& os, const SparseAlleleSet& alleles) noexcept;
        bool operator==(const SparseAlleleSet& rhs) const;
        bool operator!=(const SparseAlleleSet& rhs) const;


        [[nodiscard]] std::string serialize() const noexcept;
        [[nodiscard]] constexpr unsigned int totalPositiveCount() const noexcept {
            return set_count_;
        }
        [[nodiscard]] constexpr unsigned int totalNegativeCount() const noexcept {
            return total_alleles_ - set_count_;
        }

        constexpr static unsigned int truePositiveCount(const SparseAlleleSet& parent, const SparseAlleleSet& child) noexcept {
            unsigned int count = 0;
            for (unsigned int i = 0; i < parent.total_alleles_; ++i) {
                if (parent.allele(i) && child.allele(i)) {
                    ++count;
                }
            }
            return count;
        }

        constexpr static unsigned int trueNegativeCount(const SparseAlleleSet& parent, const SparseAlleleSet& child) noexcept {
            unsigned int count = 0;
            for (unsigned int i = 0; i < parent.total_alleles_; ++i) {
                if (!parent.allele(i) && !child.allele(i)) {
                    ++count;
                }
            }
            return count;
        }

        constexpr static unsigned int falsePositiveCount(const SparseAlleleSet& parent, const SparseAlleleSet& child) noexcept {
            unsigned int count = 0;
            for (unsigned int i = 0; i < parent.total_alleles_; ++i) {
                if (!parent.allele(i) && child.allele(i)) {
                    ++count;
                }
            }
            return count;
        }

        constexpr static unsigned int falseNegativeCount(const SparseAlleleSet& parent, const SparseAlleleSet& child) noexcept {
            unsigned int count = 0;
            for (unsigned int i = 0; i < parent.total_alleles_; ++i) {
                if (parent.allele(i) && !child.allele(i)) {
                    ++count;
                }
            }
            return count;
        }

        constexpr static SparseAlleleSet shared(const SparseAlleleSet& lhs, const SparseAlleleSet& rhs) noexcept {
            unsigned char total_alleles_ = std::min(lhs.total_alleles_, rhs.total_alleles_);
            unsigned char set_count_ = 0;
            std::array<unsigned char, MAX_SIZE> set_positions_{};
            for (unsigned int i = 0; i < total_alleles_; ++i) {
                if (lhs.allele(i) && rhs.allele(i)) {
                    set_positions_[set_count_] = i;
                    ++set_count_;
                }
            }
            return SparseAlleleSet{total_alleles_, set_positions_};
        }

        constexpr static SparseAlleleSet any(const SparseAlleleSet& lhs, const SparseAlleleSet& rhs) noexcept {
            unsigned char total_alleles_ = std::max(lhs.total_alleles_, rhs.total_alleles_);
            unsigned char set_count_ = 0;
            std::array<unsigned char, MAX_SIZE> set_positions_{};
            for (unsigned int i = 0; i < lhs.total_alleles_; ++i) {
                if (lhs.allele(i) || rhs.allele(i)) {
                    set_positions_[set_count_] = i;
                    ++set_count_;
                }
            }
            return SparseAlleleSet{total_alleles_, set_positions_};
        }

        constexpr static SparseAlleleSet diff(const SparseAlleleSet& lhs, const SparseAlleleSet& rhs) noexcept {
            unsigned char total_alleles_ = lhs.total_alleles_;
            unsigned char set_count_ = 0;
            std::array<unsigned char, MAX_SIZE> set_positions_{};
            for (unsigned int i = 0; i < lhs.total_alleles_; ++i) {
                if ((lhs.allele(i) &&!rhs.allele(i)) || (!lhs.allele(i) && rhs.allele(i))) {
                    set_positions_[set_count_] = i;
                    ++set_count_;
                }
            }
            return SparseAlleleSet{total_alleles_, set_positions_};
        }

        constexpr static SparseAlleleSet invert(const SparseAlleleSet& alleles) noexcept {
            unsigned char total_alleles_ = alleles.total_alleles_;
            unsigned char set_count_ = 0;
            std::array<unsigned char, MAX_SIZE> set_positions_{};
            for (unsigned int i = 0; i < alleles.total_alleles_; ++i) {
                if (!alleles.allele(i)) {
                    set_positions_[set_count_] = i;
                    ++set_count_;
                }
            }
            return SparseAlleleSet{total_alleles_, set_positions_};
        }

        [[nodiscard]] std::string allelesStr() const noexcept;
        [[nodiscard]] std::string compactAllelesStr() const noexcept;

        [[nodiscard]] constexpr unsigned int totalAlleles() const noexcept {
            return total_alleles_;
        }

        constexpr void flip(const size_t pos) noexcept {
            for (unsigned int i = 0; i < set_count_; ++i) {
                if (set_positions_[i] == pos) {
                    set_positions_[i] = set_positions_[--set_count_];
                    return;
                }
            }

            set_positions_[set_count_] = pos;
            ++set_count_;
            // assert(set_count_ < MAX_SIZE);
        }
        // void flip(std::vector<unsigned int> pos) noexcept;

        constexpr void set(const size_t pos, const bool val) noexcept {
            for (unsigned int i = 0; i < set_count_; ++i) {
                if (set_positions_[i] == pos) {
                    if (!val) {
                        set_positions_[i] = set_positions_[--set_count_];
                    }
                    return;
                }
            }

            if (val) {
                set_positions_[set_count_] = pos;
                ++set_count_;
                // assert(set_count_ < MAX_SIZE);
            }
        }
        // void set() noexcept;

        constexpr void reset(const size_t pos) noexcept {
            for (unsigned int i = 0; i < set_count_; ++i) {
                if (set_positions_[i] == pos) {
                    set_positions_[i] = set_positions_[--set_count_];
                    break;
                }
            }
        }
        constexpr void reset() noexcept {
            set_count_ = 0;
        }

        [[nodiscard]] constexpr SparseAlleleSet mutationMask(const SparseAlleleSet& parent) const noexcept {
            std::array<unsigned char, MAX_SIZE> mask_positions{};
            unsigned char local_set_count_ = 0;
            for (const auto pos : set_positions_) {
                if (!parent.allele(pos)) {
                    mask_positions[local_set_count_] = pos;
                    ++local_set_count_;
                }
            }

            return SparseAlleleSet(total_alleles_, mask_positions);
        }


        [[nodiscard]] constexpr bool allele(const size_t pos) const noexcept {
            for (unsigned int i = 0; i < set_count_; ++i) {
                if (set_positions_[i] == pos) {
                    return true;
                }
            }
            return false;
        }



    private:
        unsigned char total_alleles_ = 0;
        unsigned char set_count_ = 0;
        unsigned char set_positions_[MAX_SIZE];

        constexpr void updateSetPositionsFromBitString(const std::string& bitstr) {
        set_count_ = 0;
        for (unsigned int i = 0; i < bitstr.size(); ++i) {
            if (bitstr[i] == '1') {
                set_positions_[set_count_] = i;
                ++set_count_;
            }
        }
    }
    };
}




#endif //SPARSEALLELESET_H
