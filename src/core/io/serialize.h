//
// Created by Maxwell Murphy on 5/26/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_SERIALIZE_H
#define TRANSMISSION_NETWORKS_APP_SERIALIZE_H

#include <string>

#include "core/containers/Infection.h"
#include "core/datatypes/Alleles.h"
#include "core/datatypes/Simplex.h"
#include "core/parameters/TransmissionNetwork.h"
#include "core/computation/OrderDerivedParentSet.h"


namespace transmission_nets::core::io {

    inline std::string serialize(const double val) noexcept {
        return fmt::format("{}", val);
    }

    inline std::string serialize(const long double val) noexcept {
        return fmt::format("{}", val);
    }

    inline std::string serialize(const long int val) noexcept {
        return fmt::format("{}", val);
    }

    inline std::string serialize(const int &val) noexcept {
        return fmt::format("{}", val);
    }

    inline std::string serialize(const unsigned int val) noexcept {
        return fmt::format("{}", val);
    }

    inline std::string serialize(const long long val) noexcept {
        return fmt::format("{}", val);
    }

    inline std::string serialize(const std::string& val) noexcept {
        return val;
    }

    template<int MAX_COI>
    std::string serialize(const datatypes::AllelesBitSet<MAX_COI>& val) noexcept {
        return val.serialize();
    }

    template<int MAX_COI>
    std::string serialize(datatypes::AllelesBitSet<MAX_COI>& val) noexcept {
        return val.serialize();
    }

    inline std::string serialize(const datatypes::Simplex& val) noexcept {
        return val.serialize();
    }

    template<typename GeneticImpl, typename LocusImpl = containers::Locus>
    std::string serialize(const std::shared_ptr<containers::Infection<GeneticImpl, LocusImpl>> val) {
        return val->serialize();
    }

    template<typename T>
    std::string serialize(const containers::ParentSet<T> ps) {
        std::string out = "{";
        for (const auto& p : ps) {
            out += serialize(p);
            out += ";";
        }
        out.pop_back();
        out += "}";
        return out;
    }


    template<typename T, typename U>
    std::string serialize(computation::OrderDerivedParentSet<T, U> &ps) {
        std::string out = "{";
        for (const auto& p : ps.value()) {
            out += serialize(*p);
            out += ";";
        }
        out.pop_back();
        out += "}";
        return out;
    }

    template<typename T>
    std::string serialize(const boost::container::flat_set<T> ps) {
        std::string out = "{";
        for (const auto& p : ps) {
            out += serialize(p);
            out += ";";
        }
        out.pop_back();
        out += "}";
        return out;
    }


    template<typename T>
    std::string serialize(const std::vector<std::shared_ptr<T>> val) {
        std::string out;
        for (const auto& el : val) {
            out += serialize(el);
            out += ",";
        }
        out.pop_back();
        return out;
    }

    template<typename T>
    std::string serialize(const std::vector<T> val) {
        std::string out;
        for (const auto& el : val) {
            out += serialize(el);
            out += ",";
        }
        out.pop_back();
        return out;
    }

    template<typename T, long unsigned int SIZE>
    std::string serialize(std::array<T, SIZE> &val) {
        std::string out;
        for (const auto& el : val) {
            out += serialize(el);
            out += ",";
        }
        out.pop_back();
        return out;
    }
}// namespace transmission_nets::core::io


#endif//TRANSMISSION_NETWORKS_APP_SERIALIZE_H
