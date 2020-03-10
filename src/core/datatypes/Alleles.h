//
// Created by Maxwell Murphy on 12/31/19.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELE_H
#define TRANSMISSION_NETWORKS_APP_ALLELE_H

#include <bitset>
#include <ostream>

template<int MaxAlleles>
class AllelesBitSet {
public:
    explicit AllelesBitSet(const std::string &bitstr) : total_alleles_(bitstr.size()), alleles_(bitstr) {
        assert(bitstr.size() <= MaxAlleles);
    }

    AllelesBitSet(const AllelesBitSet<MaxAlleles>& other) : total_alleles_(other.total_alleles_), alleles_(other.alleles_) {
        std::cout << "alleles copy c'tor" << std::endl;
    };

    friend std::ostream &operator<<(std::ostream &os, const AllelesBitSet &alleles) noexcept {
        os << alleles.allelesStr();
        return os;
    };

    [[nodiscard]] constexpr inline unsigned int
    totalPositiveCount() noexcept {
        return alleles_.count();
    }

    [[nodiscard]] constexpr inline static unsigned int
    truePositiveCount(const AllelesBitSet<MaxAlleles> &parent, const AllelesBitSet<MaxAlleles> &child) noexcept {
        return (child.alleles_ & parent.alleles_).count();
    }

    [[nodiscard]] constexpr inline static unsigned int
    trueNegativeCount(const AllelesBitSet<MaxAlleles> &parent, const AllelesBitSet<MaxAlleles> &child) noexcept {
        assert(child.totalAlleles() == parent.totalAlleles());
        return (~child.alleles_ & ~parent.alleles_).count() - (MaxAlleles - child.totalAlleles());
    };

    [[nodiscard]] constexpr inline static unsigned int
    falsePositiveCount(const AllelesBitSet<MaxAlleles> &parent, const AllelesBitSet<MaxAlleles> &child) noexcept {
        return (child.alleles_ & ~parent.alleles_).count();
    };

    [[nodiscard]] constexpr inline static unsigned int
    falseNegativeCount(const AllelesBitSet<MaxAlleles> &parent, const AllelesBitSet<MaxAlleles> &child) noexcept {
        return (~child.alleles_ & parent.alleles_).count();
    }

    [[nodiscard]] inline std::string allelesStr() const noexcept {
        return alleles_.to_string().substr(MaxAlleles - total_alleles_, total_alleles_);
    };

    [[nodiscard]] constexpr int totalAlleles() const noexcept {
        return total_alleles_;
    };

    inline constexpr void flip(size_t pos) noexcept {
        assert(pos < total_alleles_);
        alleles_.flip(total_alleles_ - 1 - pos);
    };

    inline constexpr void set(size_t pos) noexcept {
        assert(pos < total_alleles_);
        alleles_.set(total_alleles_ - 1 - pos);
    };

    inline constexpr void reset(size_t pos) noexcept {
        assert(pos < total_alleles_);
        alleles_.reset(total_alleles_ - 1 - pos);
    };

    inline constexpr bool allele(size_t pos) noexcept {
        // Bitsets are accessed right to left so we're converting to left to right accession
        assert(pos < total_alleles_);
        return alleles_[total_alleles_ - 1 - pos];
    };

private:
    unsigned int total_alleles_;
    std::bitset<MaxAlleles> alleles_;
};


#endif //TRANSMISSION_NETWORKS_APP_ALLELE_H