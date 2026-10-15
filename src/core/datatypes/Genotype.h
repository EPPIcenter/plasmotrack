//
// A Genotype: the set of Alleles at one Locus, as a value. What an Observed genotype and a Latent genotype each hold.
// Alleles are numbered 0 to totalAlleles() - 1 in the Locus's order; a string of '0' and '1', one per Allele, reads
// and writes one. The set algebra (shared, any, diff, invert, the pairwise counts) is over the Locus's Alleles, and
// the value hashes and compares by its Alleles. Holds at most MAX_ALLELES Alleles, the engine's one limit.
//

#ifndef TRANSMISSION_NETWORKS_APP_GENOTYPE_H
#define TRANSMISSION_NETWORKS_APP_GENOTYPE_H

#include <bitset>
#include <cassert>
#include <cstddef>
#include <functional>
#include <ostream>
#include <string>
#include <vector>


namespace transmission_nets::core::datatypes {
    class Genotype {
    public:
        static constexpr int MAX_ALLELES = 128;

        explicit Genotype(const std::string& bitstr);
        explicit Genotype(int totalAlleles);
        Genotype();

        friend std::ostream& operator<<(std::ostream& os, const Genotype& alleles) noexcept;
        bool operator==(const Genotype& rhs) const;
        bool operator!=(const Genotype& rhs) const;

        Genotype(const Genotype& other)
            : total_alleles_(other.total_alleles_),
              alleles_(other.alleles_),
              mask_(other.mask_) {
            // fmt::print("Copy constructor called\n");
        }

        Genotype(Genotype&& other) noexcept
            : total_alleles_(other.total_alleles_),
              alleles_(std::move(other.alleles_)),
              mask_(std::move(other.mask_)) {
            // fmt::print("Move constructor called\n");
        }

        Genotype& operator=(const Genotype& other) {
            // fmt::print("Copy assignment called\n");
            if (this == &other)
                return *this;
            total_alleles_ = other.total_alleles_;
            alleles_       = other.alleles_;
            mask_          = other.mask_;
            return *this;
        }

        Genotype& operator=(Genotype&& other) noexcept {
            // fmt::print("Move assignment called\n");
            if (this == &other)
                return *this;
            total_alleles_ = other.total_alleles_;
            alleles_       = std::move(other.alleles_);
            mask_          = std::move(other.mask_);
            return *this;
        }

        [[nodiscard]] std::string serialize() const noexcept;

        [[nodiscard]] constexpr unsigned int
        totalPositiveCount() const noexcept;

        [[nodiscard]] constexpr static unsigned int
        truePositiveCount(const Genotype& parent, const Genotype& child) noexcept;

        [[nodiscard]] constexpr static unsigned int
        trueNegativeCount(const Genotype& parent, const Genotype& child) noexcept;

        [[nodiscard]] constexpr static unsigned int
        falsePositiveCount(const Genotype& parent, const Genotype& child) noexcept;

        [[nodiscard]] constexpr static unsigned int
        falseNegativeCount(const Genotype& parent, const Genotype& child) noexcept;

        [[nodiscard]] static Genotype
                shared(const Genotype& lhs, const Genotype& rhs) noexcept;

        [[nodiscard]] static Genotype
                any(const Genotype& lhs, const Genotype& rhs) noexcept;

        [[nodiscard]] static Genotype
                diff(const Genotype& lhs, const Genotype& rhs) noexcept;

        [[nodiscard]] static Genotype
                invert(const Genotype& alleles) noexcept;

        [[nodiscard]] std::string allelesStr() const noexcept;

        [[nodiscard]] constexpr unsigned int totalAlleles() const noexcept;

        constexpr void flip(size_t pos) noexcept;
        void flip(std::vector<unsigned int> pos) noexcept;

        constexpr void set(size_t pos, bool val = true) noexcept;
        constexpr void set() noexcept;

        constexpr void reset(size_t pos) noexcept;
        constexpr void reset() noexcept;

        Genotype mutationMask(const Genotype& parent) const noexcept;

        [[nodiscard]] constexpr bool allele(size_t pos) const noexcept;

        // Consistent with operator==, which compares total_alleles_ and alleles_.
        [[nodiscard]] std::size_t hash() const noexcept {
            return std::hash<std::bitset<MAX_ALLELES>>{}(alleles_) ^ (static_cast<std::size_t>(total_alleles_) * 0x9e3779b97f4a7c15ULL);
        }

    private:
        unsigned int total_alleles_ = 0;
        std::bitset<MAX_ALLELES> alleles_{};
        std::bitset<MAX_ALLELES> mask_{};
    };

    inline Genotype::Genotype(const std::string& bitstr) : total_alleles_(bitstr.size()), alleles_(bitstr) {
        for (size_t i = 0; i < total_alleles_; ++i) {
            mask_.set(i);
        }
        assert(bitstr.size() <= MAX_ALLELES);
    }


    inline Genotype::Genotype(int totalAlleles) : total_alleles_(totalAlleles) {
        // As the string of totalAlleles - 1 zeros and a one: the last allele present, which is bit 0.
        assert(total_alleles_ <= MAX_ALLELES);
        alleles_.reset();
        alleles_.set(0);
    }

    inline std::ostream& operator<<(std::ostream& os, const Genotype& alleles) noexcept {
        os << alleles.allelesStr();
        return os;
    }

    constexpr unsigned int Genotype::totalPositiveCount() const noexcept {
        return alleles_.count();
    }


    constexpr unsigned int Genotype::truePositiveCount(const Genotype& parent,
                                                                        const Genotype& child) noexcept {
        return (child.alleles_ & parent.alleles_).count();
    }

    constexpr unsigned int Genotype::trueNegativeCount(const Genotype& parent,
                                                                        const Genotype& child) noexcept {
        return (~child.alleles_ & ~parent.alleles_).count() - (MAX_ALLELES - child.totalAlleles());
        //        return ~(child.alleles_ | parent.alleles_).count() - (MAX_ALLELES - child.totalAlleles());
    }

    constexpr unsigned int Genotype::falsePositiveCount(const Genotype& parent,
                                                                         const Genotype& child) noexcept {
        return (child.alleles_ & ~parent.alleles_).count();
    }

    constexpr unsigned int Genotype::falseNegativeCount(const Genotype& parent,
                                                                         const Genotype& child) noexcept {

        return (~child.alleles_ & parent.alleles_).count();
    }

    inline std::string Genotype::allelesStr() const noexcept {
        return alleles_.to_string().substr(MAX_ALLELES - total_alleles_, total_alleles_);
    }

    constexpr unsigned int Genotype::totalAlleles() const noexcept {
        return total_alleles_;
    }

    constexpr void Genotype::flip(size_t pos) noexcept {
        assert(pos < total_alleles_);
        alleles_.flip(total_alleles_ - 1 - pos);
    }

    inline void Genotype::flip(std::vector<unsigned int> pos) noexcept {
        for (auto& p : pos) {
            flip(p);
        }
    }

    constexpr void Genotype::set(size_t pos, bool val) noexcept {
        assert(pos < total_alleles_);
        alleles_.set(total_alleles_ - 1 - pos, val);
    }

    constexpr void Genotype::set() noexcept {
        alleles_.set();
    }

    constexpr void Genotype::reset(size_t pos) noexcept {
        assert(pos < total_alleles_);
        alleles_.reset(total_alleles_ - 1 - pos);
    }

    constexpr void Genotype::reset() noexcept {
        alleles_.reset();
    }

    constexpr bool Genotype::allele(size_t pos) const noexcept {
        // Bitsets are accessed right to left so we're converting to left to right accession
        assert(pos < total_alleles_);
        assert((total_alleles_ - 1 - pos) >= 0);
        return alleles_[total_alleles_ - 1 - pos];
    }

    inline std::string Genotype::serialize() const noexcept {
        return allelesStr();
    }

    inline Genotype::Genotype() {
        total_alleles_ = 0;
        alleles_       = std::bitset<MAX_ALLELES>{};
    }

    inline Genotype Genotype::shared(const Genotype& lhs, const Genotype& rhs) noexcept {
        Genotype result;
        result.alleles_ = lhs.alleles_ & rhs.alleles_;
        result.total_alleles_ = lhs.total_alleles_;
        result.mask_ = lhs.mask_;
        return result;
    }

    inline Genotype Genotype::any(const Genotype& lhs, const Genotype& rhs) noexcept {
        Genotype result;
        result.alleles_ = lhs.alleles_ | rhs.alleles_;
        result.total_alleles_ = lhs.total_alleles_;
        result.mask_ = lhs.mask_;
        return result;
    }

    inline Genotype Genotype::diff(const Genotype& lhs, const Genotype& rhs) noexcept {
        Genotype result;
        result.alleles_ = lhs.alleles_ ^ rhs.alleles_;
        result.total_alleles_ = lhs.total_alleles_;
        result.mask_ = lhs.mask_;
        return result;
    }

    inline Genotype Genotype::invert(const Genotype& alleles) noexcept {
        Genotype result;
        result.alleles_ = alleles.alleles_;
        result.total_alleles_ = alleles.total_alleles_;
        result.mask_ = alleles.mask_;
        result.alleles_ = ~result.alleles_ & result.mask_;
        return result;
    }



    inline bool Genotype::operator==(const Genotype& rhs) const {
        return total_alleles_ == rhs.total_alleles_ &&
               alleles_ == rhs.alleles_;
    }

    inline bool Genotype::operator!=(const Genotype& rhs) const {
        return !(*this == rhs);
    }

    inline Genotype Genotype::mutationMask(const Genotype& parent) const noexcept {
        // The mutation mask gets the alleles that are present in the child but not in the parent
        Genotype result;
        result.alleles_ = ~parent.alleles_ & this->alleles_;
        result.total_alleles_ = this->total_alleles_;
        result.mask_ = this->mask_;
        return result;
    }
}// namespace transmission_nets::core::datatypes


#endif//TRANSMISSION_NETWORKS_APP_GENOTYPE_H