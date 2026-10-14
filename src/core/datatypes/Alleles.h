//
// Created by Maxwell Murphy on 12/31/19.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELE_H
#define TRANSMISSION_NETWORKS_APP_ALLELE_H

#include <bitset>
#include <cassert>
#include <vector>


namespace transmission_nets::core::datatypes {
    template<int MaxAlleles>
    class AllelesBitSet {
    public:
        explicit AllelesBitSet(const std::string& bitstr);
        explicit AllelesBitSet(int totalAlleles);
        AllelesBitSet();

        template<int T>
        friend std::ostream& operator<<(std::ostream& os, const AllelesBitSet<T>& alleles) noexcept;
        bool operator==(const AllelesBitSet& rhs) const;
        bool operator!=(const AllelesBitSet& rhs) const;

        AllelesBitSet(const AllelesBitSet& other)
            : total_alleles_(other.total_alleles_),
              alleles_(other.alleles_),
              mask_(other.mask_) {
            // fmt::print("Copy constructor called\n");
        }

        AllelesBitSet(AllelesBitSet&& other) noexcept
            : total_alleles_(other.total_alleles_),
              alleles_(std::move(other.alleles_)),
              mask_(std::move(other.mask_)) {
            // fmt::print("Move constructor called\n");
        }

        AllelesBitSet& operator=(const AllelesBitSet& other) {
            // fmt::print("Copy assignment called\n");
            if (this == &other)
                return *this;
            total_alleles_ = other.total_alleles_;
            alleles_       = other.alleles_;
            mask_          = other.mask_;
            return *this;
        }

        AllelesBitSet& operator=(AllelesBitSet&& other) noexcept {
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

        [[nodiscard]] constexpr unsigned int
        totalNegativeCount() const noexcept;

        [[nodiscard]] constexpr static unsigned int
        truePositiveCount(const AllelesBitSet<MaxAlleles>& parent, const AllelesBitSet<MaxAlleles>& child) noexcept;

        [[nodiscard]] constexpr static unsigned int
        trueNegativeCount(const AllelesBitSet<MaxAlleles>& parent, const AllelesBitSet<MaxAlleles>& child) noexcept;

        [[nodiscard]] constexpr static unsigned int
        falsePositiveCount(const AllelesBitSet<MaxAlleles>& parent, const AllelesBitSet<MaxAlleles>& child) noexcept;

        [[nodiscard]] constexpr static unsigned int
        falseNegativeCount(const AllelesBitSet<MaxAlleles>& parent, const AllelesBitSet<MaxAlleles>& child) noexcept;

        [[nodiscard]] constexpr static AllelesBitSet<MaxAlleles>
                shared(const AllelesBitSet<MaxAlleles>& lhs, const AllelesBitSet<MaxAlleles>& rhs) noexcept;

        [[nodiscard]] constexpr static AllelesBitSet<MaxAlleles>
                any(const AllelesBitSet<MaxAlleles>& lhs, const AllelesBitSet<MaxAlleles>& rhs) noexcept;

        [[nodiscard]] constexpr static AllelesBitSet<MaxAlleles>
                diff(const AllelesBitSet<MaxAlleles>& lhs, const AllelesBitSet<MaxAlleles>& rhs) noexcept;

        [[nodiscard]] constexpr static AllelesBitSet<MaxAlleles>
                invert(const AllelesBitSet<MaxAlleles>& alleles) noexcept;

        [[nodiscard]] std::string allelesStr() const noexcept;
        [[nodiscard]] std::string compactAllelesStr() const noexcept;

        [[nodiscard]] constexpr unsigned int totalAlleles() const noexcept;

        constexpr void flip(size_t pos) noexcept;
        void flip(std::vector<unsigned int> pos) noexcept;

        constexpr void set(size_t pos, bool val = true) noexcept;
        constexpr void set() noexcept;

        constexpr void reset(size_t pos) noexcept;
        constexpr void reset() noexcept;

        constexpr AllelesBitSet<MaxAlleles> mutationMask(const AllelesBitSet<MaxAlleles>& parent) const noexcept;

        [[nodiscard]] constexpr bool allele(size_t pos) const noexcept;

        // Consistent with operator==, which compares total_alleles_ and alleles_.
        [[nodiscard]] std::size_t hash() const noexcept {
            return std::hash<std::bitset<MaxAlleles>>{}(alleles_) ^ (static_cast<std::size_t>(total_alleles_) * 0x9e3779b97f4a7c15ULL);
        }

    private:
        unsigned int total_alleles_ = 0;
        std::bitset<MaxAlleles> alleles_{};
        std::bitset<MaxAlleles> mask_{};
    };

    template<int MaxAlleles>
    AllelesBitSet<MaxAlleles>::AllelesBitSet(const std::string& bitstr) : total_alleles_(bitstr.size()), alleles_(bitstr) {
        for (size_t i = 0; i < total_alleles_; ++i) {
            mask_.set(i);
        }
        assert(bitstr.size() <= MaxAlleles);
    }


    template<int MaxAlleles>
    AllelesBitSet<MaxAlleles>::AllelesBitSet(int totalAlleles) : total_alleles_(totalAlleles) {
        // As the string of totalAlleles - 1 zeros and a one: the last allele present, which is bit 0.
        assert(total_alleles_ <= MaxAlleles);
        alleles_.reset();
        alleles_.set(0);
    }

    template<int MaxAlleles>
    std::ostream& operator<<(std::ostream& os, const AllelesBitSet<MaxAlleles>& alleles) noexcept {
        os << alleles.allelesStr();
        return os;
    }

    template<int MaxAlleles>
    constexpr unsigned int AllelesBitSet<MaxAlleles>::totalPositiveCount() const noexcept {
        return alleles_.count();
    }

    template<int MaxAlleles>
    constexpr unsigned int AllelesBitSet<MaxAlleles>::totalNegativeCount() const noexcept {
        return total_alleles_ - alleles_.count();
    }


    template<int MaxAlleles>
    constexpr unsigned int AllelesBitSet<MaxAlleles>::truePositiveCount(const AllelesBitSet<MaxAlleles>& parent,
                                                                        const AllelesBitSet<MaxAlleles>& child) noexcept {
        return (child.alleles_ & parent.alleles_).count();
    }

    template<int MaxAlleles>
    constexpr unsigned int AllelesBitSet<MaxAlleles>::trueNegativeCount(const AllelesBitSet<MaxAlleles>& parent,
                                                                        const AllelesBitSet<MaxAlleles>& child) noexcept {
        return (~child.alleles_ & ~parent.alleles_).count() - (MaxAlleles - child.totalAlleles());
        //        return ~(child.alleles_ | parent.alleles_).count() - (MaxAlleles - child.totalAlleles());
    }

    template<int MaxAlleles>
    constexpr unsigned int AllelesBitSet<MaxAlleles>::falsePositiveCount(const AllelesBitSet<MaxAlleles>& parent,
                                                                         const AllelesBitSet<MaxAlleles>& child) noexcept {
        return (child.alleles_ & ~parent.alleles_).count();
    }

    template<int MaxAlleles>
    constexpr unsigned int AllelesBitSet<MaxAlleles>::falseNegativeCount(const AllelesBitSet<MaxAlleles>& parent,
                                                                         const AllelesBitSet<MaxAlleles>& child) noexcept {

        return (~child.alleles_ & parent.alleles_).count();
    }

    template<int MaxAlleles>
    std::string AllelesBitSet<MaxAlleles>::allelesStr() const noexcept {
        return alleles_.to_string().substr(MaxAlleles - total_alleles_, total_alleles_);
    }

    template<int MaxAlleles>
    std::string AllelesBitSet<MaxAlleles>::compactAllelesStr() const noexcept {
        return std::to_string(std::stoll(allelesStr()));
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
    void AllelesBitSet<MaxAlleles>::flip(std::vector<unsigned int> pos) noexcept {
        for (auto& p : pos) {
            flip(p);
        }
    }

    template<int MaxAlleles>
    constexpr void AllelesBitSet<MaxAlleles>::set(size_t pos, bool val) noexcept {
        assert(pos < total_alleles_);
        alleles_.set(total_alleles_ - 1 - pos, val);
    }

    template<int MaxAlleles>
    constexpr void AllelesBitSet<MaxAlleles>::set() noexcept {
        alleles_.set();
    }

    template<int MaxAlleles>
    constexpr void AllelesBitSet<MaxAlleles>::reset(size_t pos) noexcept {
        assert(pos < total_alleles_);
        alleles_.reset(total_alleles_ - 1 - pos);
    }

    template<int MaxAlleles>
    constexpr void AllelesBitSet<MaxAlleles>::reset() noexcept {
        alleles_.reset();
    }

    template<int MaxAlleles>
    constexpr bool AllelesBitSet<MaxAlleles>::allele(size_t pos) const noexcept {
        // Bitsets are accessed right to left so we're converting to left to right accession
        assert(pos < total_alleles_);
        assert((total_alleles_ - 1 - pos) >= 0);
        return alleles_[total_alleles_ - 1 - pos];
    }

    template<int MaxAlleles>
    std::string AllelesBitSet<MaxAlleles>::serialize() const noexcept {
        return allelesStr();
    }

    template<int MaxAlleles>
    AllelesBitSet<MaxAlleles>::AllelesBitSet() {
        total_alleles_ = 0;
        alleles_       = std::bitset<MaxAlleles>{};
    }

    template<int MaxAlleles>
    constexpr AllelesBitSet<MaxAlleles> AllelesBitSet<MaxAlleles>::shared(const AllelesBitSet<MaxAlleles>& lhs, const AllelesBitSet<MaxAlleles>& rhs) noexcept {
        AllelesBitSet<MaxAlleles> result;
        result.alleles_ = lhs.alleles_ & rhs.alleles_;
        result.total_alleles_ = lhs.total_alleles_;
        result.mask_ = lhs.mask_;
        return result;
    }

    template<int MaxAlleles>
    constexpr AllelesBitSet<MaxAlleles> AllelesBitSet<MaxAlleles>::any(const AllelesBitSet<MaxAlleles>& lhs, const AllelesBitSet<MaxAlleles>& rhs) noexcept {
        AllelesBitSet<MaxAlleles> result;
        result.alleles_ = lhs.alleles_ | rhs.alleles_;
        result.total_alleles_ = lhs.total_alleles_;
        result.mask_ = lhs.mask_;
        return result;
    }

    template<int MaxAlleles>
    constexpr AllelesBitSet<MaxAlleles> AllelesBitSet<MaxAlleles>::diff(const AllelesBitSet<MaxAlleles>& lhs, const AllelesBitSet<MaxAlleles>& rhs) noexcept {
        AllelesBitSet<MaxAlleles> result;
        result.alleles_ = lhs.alleles_ ^ rhs.alleles_;
        result.total_alleles_ = lhs.total_alleles_;
        result.mask_ = lhs.mask_;
        return result;
    }

    template<int MaxAlleles>
    constexpr AllelesBitSet<MaxAlleles> AllelesBitSet<MaxAlleles>::invert(const AllelesBitSet<MaxAlleles>& alleles) noexcept {
        AllelesBitSet<MaxAlleles> result;
        result.alleles_ = alleles.alleles_;
        result.total_alleles_ = alleles.total_alleles_;
        result.mask_ = alleles.mask_;
        result.alleles_ = ~result.alleles_ & result.mask_;
        return result;
    }



    template<int MaxAlleles>
    bool AllelesBitSet<MaxAlleles>::operator==(const AllelesBitSet& rhs) const {
        return total_alleles_ == rhs.total_alleles_ &&
               alleles_ == rhs.alleles_;
    }

    template<int MaxAlleles>
    bool AllelesBitSet<MaxAlleles>::operator!=(const AllelesBitSet& rhs) const {
        return !(*this == rhs);
    }

    template<int MaxAlleles>
    constexpr AllelesBitSet<MaxAlleles> AllelesBitSet<MaxAlleles>::mutationMask(const AllelesBitSet<MaxAlleles>& parent) const noexcept {
        // The mutation mask gets the alleles that are present in the child but not in the parent
        AllelesBitSet<MaxAlleles> result;
        result.alleles_ = ~parent.alleles_ & this->alleles_;
        result.total_alleles_ = this->total_alleles_;
        result.mask_ = this->mask_;
        return result;
    }
}// namespace transmission_nets::core::datatypes


#endif//TRANSMISSION_NETWORKS_APP_ALLELE_H