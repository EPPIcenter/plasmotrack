//
// Created by Maxwell Murphy on 12/30/19.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELES_H
#define TRANSMISSION_NETWORKS_APP_ALLELES_H

template <int MaxAlleles>
class Alleles : public Parameter<std::bitset<MaxAlleles>> {
private:
    int num_alleles_{};
    Alleles(const std::string &id, const std::bitset<MaxAlleles> &value);

public:
    Alleles(const std::string &id, const std::string &bitstr);

    int num_alleles();

};

template<int MaxAlleles>
Alleles<MaxAlleles>::Alleles(const std::string &id, const std::bitset<MaxAlleles> &value) : Parameter<std::bitset<MaxAlleles>>(id, value) {}

template<int MaxAlleles>
Alleles<MaxAlleles>::Alleles(const std::string &id, const std::string &bitstr):Parameter<std::bitset<MaxAlleles>>(id) {
    assert(bitstr.size() <= MaxAlleles);
    this->value_ = std::bitset<MaxAlleles>(bitstr);
    this->num_alleles_ = bitstr.size();
    this->setClean();
}

template<int MaxAlleles>
int Alleles<MaxAlleles>::num_alleles() {
    return this->num_alleles_;
}


#endif //TRANSMISSION_NETWORKS_APP_ALLELES_H
