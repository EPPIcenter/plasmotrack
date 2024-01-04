//
// Created by Maxwell Murphy on 5/25/23.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELESBITVEC_H
#define TRANSMISSION_NETWORKS_APP_ALLELESBITVEC_H

#include <numeric>
#include <string>
#include <vector>

namespace transmission_nets::core::datatypes {
    std::vector<bool> bitstr_to_bitvec(const std::string& bitstr) {
        std::vector<bool> bitvec;
        bitvec.reserve(bitstr.size());
        for (const char& c : bitstr) {
            bitvec.push_back(c == '1');
        }
        return bitvec;
    };

    class AllelesBitVec {
    public:
        explicit AllelesBitVec(const std::string& bitstr);
        AllelesBitVec(int totalAlleles);
        AllelesBitVec();

        template<int T>
        friend std::ostream& operator<<(std::ostream& os, const AllelesBitVec& alleles) noexcept;
        bool operator==(const AllelesBitVec& rhs) const;
        bool operator!=(const AllelesBitVec& rhs) const;

        [[nodiscard]] std::string serialize() const noexcept;

        [[nodiscard]] constexpr inline unsigned int
        totalPositiveCount() const noexcept;

        [[nodiscard]] constexpr inline unsigned int
        totalNegativeCount() const noexcept;

        [[nodiscard]] constexpr inline static unsigned int
        truePositiveCount(const AllelesBitVec& parent, const AllelesBitVec& child) noexcept;

        [[nodiscard]] constexpr inline static unsigned int
        trueNegativeCount(const AllelesBitVec& parent, const AllelesBitVec& child) noexcept;

        [[nodiscard]] constexpr inline static unsigned int
        falsePositiveCount(const AllelesBitVec& parent, const AllelesBitVec& child) noexcept;

        [[nodiscard]] constexpr inline static unsigned int
        falseNegativeCount(const AllelesBitVec& parent, const AllelesBitVec& child) noexcept;

        [[nodiscard]] constexpr inline static AllelesBitVec
        shared(const AllelesBitVec& lhs, const AllelesBitVec& rhs) noexcept;

        [[nodiscard]] constexpr inline static AllelesBitVec
        any(const AllelesBitVec& lhs, const AllelesBitVec& rhs) noexcept;

        [[nodiscard]] constexpr inline static AllelesBitVec
        diff(const AllelesBitVec& lhs, const AllelesBitVec& rhs) noexcept;

        [[nodiscard]] constexpr inline static AllelesBitVec
        invert(const AllelesBitVec& alleles) noexcept;

        [[nodiscard]] inline std::string allelesStr() const noexcept;
        [[nodiscard]] inline std::string compactAllelesStr() const noexcept;

        [[nodiscard]] constexpr unsigned int totalAlleles() const noexcept;

        inline constexpr void flip(size_t pos) noexcept;
        void flip(std::vector<unsigned int> pos) noexcept;

        inline constexpr void set(size_t pos, bool val = true) noexcept;
        inline constexpr void set() noexcept;

        inline constexpr void reset(size_t pos) noexcept;
        inline constexpr void reset() noexcept;

        inline constexpr AllelesBitVec mutationMask(const AllelesBitVec& parent) const noexcept;

        [[nodiscard]] inline constexpr bool allele(size_t pos) const noexcept;

    private:
        unsigned int total_alleles_ = 0;
        std::vector<bool> alleles_{};
        std::vector<bool> mask_{};
    };


    AllelesBitVec::AllelesBitVec(const std::string& bitstr) : total_alleles_(bitstr.size()) {
        alleles_ = bitstr_to_bitvec(bitstr);
        for (size_t i = 0; i < total_alleles_; ++i) {
            mask_[i] = true;
        }
    }



    AllelesBitVec::AllelesBitVec(int totalAlleles) : total_alleles_(totalAlleles) {
        alleles_ = std::vector(total_alleles_, false);
        alleles_[total_alleles_ - 1] = true;
    }


    std::ostream& operator<<(std::ostream& os, const AllelesBitVec& alleles) noexcept {
        os << alleles.allelesStr();
        return os;
    }


    constexpr unsigned int AllelesBitVec::totalPositiveCount() const noexcept {
        // add up all the true values in the bit vector
        return sum(alleles_.begin(), alleles_.end());
    }


    constexpr unsigned int AllelesBitVec::totalNegativeCount() const noexcept {
        return total_alleles_ - alleles_.count();
    }



    constexpr unsigned int AllelesBitVec::truePositiveCount(const AllelesBitVec& parent,
                                                                        const AllelesBitVec& child) noexcept {
        return (child.alleles_ & parent.alleles_).count();
    }


    constexpr unsigned int AllelesBitVec::trueNegativeCount(const AllelesBitVec& parent,
                                                                        const AllelesBitVec& child) noexcept {
        return (~child.alleles_ & ~parent.alleles_).count() - (MaxAlleles - child.totalAlleles());
        //        return ~(child.alleles_ | parent.alleles_).count() - (MaxAlleles - child.totalAlleles());
    }


    constexpr unsigned int AllelesBitVec::falsePositiveCount(const AllelesBitVec& parent,
                                                                         const AllelesBitVec& child) noexcept {
        return (child.alleles_ & ~parent.alleles_).count();
    }


    constexpr unsigned int AllelesBitVec::falseNegativeCount(const AllelesBitVec& parent,
                                                                         const AllelesBitVec& child) noexcept {

        return (~child.alleles_ & parent.alleles_).count();
    }


    std::string AllelesBitVec::allelesStr() const noexcept {
        return alleles_.to_string().substr(MaxAlleles - total_alleles_, total_alleles_);
    }


    std::string AllelesBitVec::compactAllelesStr() const noexcept {
        return std::to_string(std::stoll(allelesStr()));
    }


    constexpr unsigned int AllelesBitVec::totalAlleles() const noexcept {
        return total_alleles_;
    }


    constexpr void AllelesBitVec::flip(size_t pos) noexcept {
        assert(pos < total_alleles_);
        alleles_.flip(total_alleles_ - 1 - pos);
    }


    void AllelesBitVec::flip(std::vector<unsigned int> pos) noexcept {
        for (auto& p : pos) {
            flip(p);
        }
    }


    constexpr void AllelesBitVec::set(size_t pos, bool val) noexcept {
        assert(pos < total_alleles_);
        alleles_.set(total_alleles_ - 1 - pos, val);
    }


    constexpr void AllelesBitVec::set() noexcept {
        alleles_.set();
    }


    constexpr void AllelesBitVec::reset(size_t pos) noexcept {
        assert(pos < total_alleles_);
        alleles_.reset(total_alleles_ - 1 - pos);
    }


    constexpr void AllelesBitVec::reset() noexcept {
        alleles_.reset();
    }


    constexpr bool AllelesBitVec::allele(size_t pos) const noexcept {
        // Bitsets are accessed right to left so we're converting to left to right accession
        assert(pos < total_alleles_);
        return alleles_[total_alleles_ - 1 - pos];
    }


    std::string AllelesBitVec::serialize() const noexcept {
        return allelesStr();
    }

    AllelesBitVec::AllelesBitVec() {
        total_alleles_ = 0;
        alleles_       = std::bitset{};
    }


    constexpr AllelesBitVec AllelesBitVec::shared(const AllelesBitVec& lhs, const AllelesBitVec& rhs) noexcept {
        AllelesBitVec result;
        result.alleles_ = lhs.alleles_ & rhs.alleles_;
        result.total_alleles_ = lhs.total_alleles_;
        result.mask_ = lhs.mask_;
        return result;
    }


    constexpr AllelesBitVec AllelesBitVec::any(const AllelesBitVec& lhs, const AllelesBitVec& rhs) noexcept {
        AllelesBitVec result;
        result.alleles_ = lhs.alleles_ | rhs.alleles_;
        result.total_alleles_ = lhs.total_alleles_;
        result.mask_ = lhs.mask_;
        return result;
    }


    constexpr AllelesBitVec AllelesBitVec::diff(const AllelesBitVec& lhs, const AllelesBitVec& rhs) noexcept {
        AllelesBitVec result;
        result.alleles_ = lhs.alleles_ ^ rhs.alleles_;
        result.total_alleles_ = lhs.total_alleles_;
        result.mask_ = lhs.mask_;
        return result;
    }


    constexpr AllelesBitVec AllelesBitVec::invert(const AllelesBitVec& alleles) noexcept {
        AllelesBitVec result;
        result.alleles_ = alleles.alleles_;
        result.total_alleles_ = alleles.total_alleles_;
        result.mask_ = alleles.mask_;
        result.alleles_ = ~result.alleles_ & result.mask_;
        return result;
    }




    bool AllelesBitVec::operator==(const AllelesBitVec& rhs) const {
        return total_alleles_ == rhs.total_alleles_ &&
               alleles_ == rhs.alleles_;
    }


    bool AllelesBitVec::operator!=(const AllelesBitVec& rhs) const {
        return !(*this == rhs);
    }


    constexpr AllelesBitVec AllelesBitVec::mutationMask(const AllelesBitVec& parent) const noexcept {
        // The mutation mask gets the alleles that are present in the child but not in the parent
        AllelesBitVec result;
        result.alleles_ = this->alleles_;
        result.alleles_.flip();
        result.alleles_ &= parent.alleles_;
        result.alleles_ = ~parent.alleles_ & this->alleles_;
        result.total_alleles_ = this->total_alleles_;
        result.mask_ = this->mask_;
        return result;

    }
}

#endif//TRANSMISSION_NETWORKS_APP_ALLELESBITVEC_H
