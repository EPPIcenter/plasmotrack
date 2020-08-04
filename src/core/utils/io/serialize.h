//
// Created by Maxwell Murphy on 5/26/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_SERIALIZE_H
#define TRANSMISSION_NETWORKS_APP_SERIALIZE_H

#include <string>

#include "core/datatypes/Alleles.h"
#include "core/datatypes/Simplex.h"
#include "core/containers/Infection.h"
#include "core/containers/TransmissionNetwork.h"

inline std::string serialize(const double val) noexcept {
    return std::to_string(val);
}

inline std::string serialize(const int val) noexcept {
    return std::to_string(val);
}


template<int MAX_COI>
std::string serialize(const AllelesBitSet<MAX_COI>& val) noexcept {
    return val.serialize();
}

inline std::string serialize(const Simplex& val) noexcept {
    return val.serialize();
}

template<typename GeneticImpl, typename LocusImpl = Locus>
std::string serialize(const Infection<GeneticImpl, LocusImpl>* val) {
    return val->serialize();
}

template<typename NodeValueImpl>
std::string serialize(const TransmissionNetwork<NodeValueImpl>& val) {
    return val->serialize();
}


template<typename T>
std::string serialize(const std::vector<T*> val) {
    std::string out;
    for (const auto& el : val) {
        out += serialize(el);
        out += ",";
    }
    out.pop_back();
    return out;
}



#endif //TRANSMISSION_NETWORKS_APP_SERIALIZE_H
