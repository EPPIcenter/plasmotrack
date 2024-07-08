//
// Created by mmurphy on 5/28/24.
//

#include "SparseAlleleSet.h"
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

namespace transmission_nets::core::datatypes {
    SparseAlleleSet::SparseAlleleSet(const std::string& bitstr) : total_alleles_(bitstr.size()) {
        updateSetPositionsFromBitString(bitstr);
    }

    SparseAlleleSet::SparseAlleleSet(unsigned int totalAlleles) : total_alleles_(totalAlleles) {
        std::memset(set_positions_, 0, sizeof(set_positions_));
    }

    SparseAlleleSet::SparseAlleleSet(const SparseAlleleSet& other) {
        set_count_ = other.set_count_;
        total_alleles_ = other.total_alleles_;
        std::memcpy(set_positions_, other.set_positions_, sizeof(set_positions_));
    }

    SparseAlleleSet::SparseAlleleSet(SparseAlleleSet&& other) noexcept {
        set_count_ = other.set_count_;
        total_alleles_ = other.total_alleles_;
        std::memcpy(set_positions_, other.set_positions_, sizeof(set_positions_));
    }

    SparseAlleleSet::SparseAlleleSet() {
        set_count_ = 0;
        total_alleles_ = 0;
        std::memset(set_positions_, 0, sizeof(set_positions_));
    }

    SparseAlleleSet& SparseAlleleSet::operator=(const SparseAlleleSet& other) {
        if (this != &other) {
            set_count_ = other.set_count_;
            total_alleles_ = other.total_alleles_;
            std::memcpy(set_positions_, other.set_positions_, sizeof(set_positions_));
        }
        return *this;
    }

    SparseAlleleSet& SparseAlleleSet::operator=(SparseAlleleSet&& other) noexcept {
        if (this != &other) {
            set_count_ = other.set_count_;
            total_alleles_ = other.total_alleles_;
            std::memcpy(set_positions_, other.set_positions_, sizeof(set_positions_));
        }
        return *this;
    }

    std::ostream& operator<<(std::ostream& os, const SparseAlleleSet& alleles) noexcept {
        os << alleles.allelesStr();
        return os;
    }

    bool SparseAlleleSet::operator==(const SparseAlleleSet& rhs) const {
        return total_alleles_ == rhs.total_alleles_ && set_count_ == rhs.set_count_ &&
               std::memcmp(set_positions_, rhs.set_positions_, sizeof(set_positions_)) != 0;
    }

    bool SparseAlleleSet::operator!=(const SparseAlleleSet& rhs) const {
        return!(*this == rhs);
    }

    std::string SparseAlleleSet::serialize() const noexcept {
        return allelesStr();
    }

    std::string SparseAlleleSet::allelesStr() const noexcept {
        std::string str(total_alleles_, '0');
        for (unsigned int i = 0; i < set_count_; ++i) {
            str[set_positions_[i]] = '1';
        }
        return str;
    }

    std::string SparseAlleleSet::compactAllelesStr() const noexcept {
        return std::to_string(std::stoull(allelesStr(), nullptr, 2));
    }

};
