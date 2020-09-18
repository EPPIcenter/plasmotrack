//
// Created by Maxwell Murphy on 12/31/19.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELE_H
#define TRANSMISSION_NETWORKS_APP_ALLELE_H

#include <bitset>
#include <cassert>


namespace transmission_nets::core::datatypes {
    template<int MaxAlleles>
    class AllelesBitSet {
    public:
        explicit AllelesBitSet(const std::string &bitstr);

//    AllelesBitSet& operator=(AllelesBitSet other) {
//        assert(total_alleles_ == other.total_alleles_);
//        if(&other == this) {
//            return *this;
//        }
//
//        alleles_ = other.alleles_;
//        return *this;
//    }

        template<int T>
        friend std::ostream &operator<<(std::ostream &os, const AllelesBitSet<T> &alleles) noexcept;

        [[nodiscard]] std::string serialize() const noexcept;

        [[nodiscard]] constexpr inline unsigned int
        totalPositiveCount() const noexcept;

        [[nodiscard]] constexpr inline static unsigned int
        truePositiveCount(const AllelesBitSet<MaxAlleles> &parent, const AllelesBitSet<MaxAlleles> &child) noexcept;

        [[nodiscard]] constexpr inline static unsigned int
        trueNegativeCount(const AllelesBitSet<MaxAlleles> &parent, const AllelesBitSet<MaxAlleles> &child) noexcept;

        [[nodiscard]] constexpr inline static unsigned int
        falsePositiveCount(const AllelesBitSet<MaxAlleles> &parent, const AllelesBitSet<MaxAlleles> &child) noexcept;

        [[nodiscard]] constexpr inline static unsigned int
        falseNegativeCount(const AllelesBitSet<MaxAlleles> &parent, const AllelesBitSet<MaxAlleles> &child) noexcept;

        [[nodiscard]] inline std::string allelesStr() const noexcept;

        [[nodiscard]] constexpr unsigned int totalAlleles() const noexcept;

        inline constexpr void flip(size_t pos) noexcept;

        inline constexpr void set(size_t pos) noexcept;

        inline constexpr void set() noexcept;

        inline constexpr void reset(size_t pos) noexcept;

        [[nodiscard]] inline constexpr bool allele(size_t pos) const noexcept;

    private:
        unsigned int total_alleles_;
        std::bitset<MaxAlleles> alleles_;
    };

    template<int MaxAlleles>
    AllelesBitSet<MaxAlleles>::AllelesBitSet(const std::string &bitstr) : total_alleles_(bitstr.size()), alleles_(bitstr) {
        assert(bitstr.size() <= MaxAlleles);
    }

    template<int MaxAlleles>
    std::ostream &operator<<(std::ostream &os, const AllelesBitSet<MaxAlleles> &alleles) noexcept {
        os << alleles.allelesStr();
        return os;
    }

    template<int MaxAlleles>
    constexpr unsigned int AllelesBitSet<MaxAlleles>::totalPositiveCount() const noexcept {
        return alleles_.count();
    }

    template<int MaxAlleles>
    constexpr unsigned int AllelesBitSet<MaxAlleles>::truePositiveCount(const AllelesBitSet<MaxAlleles> &parent,
                                                                        const AllelesBitSet<MaxAlleles> &child) noexcept {
        return (child.alleles_ & parent.alleles_).count();
    }

    template<int MaxAlleles>
    constexpr unsigned int AllelesBitSet<MaxAlleles>::trueNegativeCount(const AllelesBitSet<MaxAlleles> &parent,
                                                                        const AllelesBitSet<MaxAlleles> &child) noexcept {
        assert(child.totalAlleles() == parent.totalAlleles());
        return (~child.alleles_ & ~parent.alleles_).count() - (MaxAlleles - child.totalAlleles());
    }

    template<int MaxAlleles>
    constexpr unsigned int AllelesBitSet<MaxAlleles>::falsePositiveCount(const AllelesBitSet<MaxAlleles> &parent,
                                                                         const AllelesBitSet<MaxAlleles> &child) noexcept {
        return (child.alleles_ & ~parent.alleles_).count();
    }

    template<int MaxAlleles>
    constexpr unsigned int AllelesBitSet<MaxAlleles>::falseNegativeCount(const AllelesBitSet<MaxAlleles> &parent,
                                                                         const AllelesBitSet<MaxAlleles> &child) noexcept {
        return (~child.alleles_ & parent.alleles_).count();
    }

    template<int MaxAlleles>
    std::string AllelesBitSet<MaxAlleles>::allelesStr() const noexcept {
        return alleles_.to_string().substr(MaxAlleles - total_alleles_, total_alleles_);
    }

    template<int MaxAlleles>
    constexpr unsigned int AllelesBitSet<MaxAlleles>::totalAlleles() const noexcept {
        return total_alleles_;
    }

    template<int MaxAlleles>
    constexpr void AllelesBitSet<MaxAlleles>::flip(size_t pos) noexcept {
        assert(pos < total_alleles_);
        alleles_.flip(total_alleles_ - 1 - pos);
    }

    template<int MaxAlleles>
    constexpr void AllelesBitSet<MaxAlleles>::set(size_t pos) noexcept {
        assert(pos < total_alleles_);
        alleles_.set(total_alleles_ - 1 - pos);
    }

    template<int MaxAlleles>
    constexpr void AllelesBitSet<MaxAlleles>::set() noexcept {
        for (int i = 0; i < total_alleles_; ++i) {
            alleles_.set(i);
        };
    }


    template<int MaxAlleles>
    constexpr void AllelesBitSet<MaxAlleles>::reset(size_t pos) noexcept {
        assert(pos < total_alleles_);
        alleles_.reset(total_alleles_ - 1 - pos);
    }

    template<int MaxAlleles>
    constexpr bool AllelesBitSet<MaxAlleles>::allele(size_t pos) const noexcept {
        // Bitsets are accessed right to left so we're converting to left to right accession
        assert(pos < total_alleles_);
        return alleles_[total_alleles_ - 1 - pos];
    }

    template<int MaxAlleles>
    std::string AllelesBitSet<MaxAlleles>::serialize() const noexcept {
        return allelesStr();
    }
}



#endif //TRANSMISSION_NETWORKS_APP_ALLELE_H