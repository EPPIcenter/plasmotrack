//
// Created by Maxwell Murphy on 2/23/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_LOCUS_H
#define TRANSMISSION_NETWORKS_APP_LOCUS_H

class Locus {
public:
    const unsigned int uid;
    const std::string label;

    explicit Locus(std::string label, int total_alleles);

    [[nodiscard]] int totalAlleles() const noexcept;

    bool operator<(const Locus &rhs) const noexcept;

    bool operator>(const Locus &rhs) const noexcept;

    bool operator<=(const Locus &rhs) const noexcept;

    bool operator>=(const Locus &rhs) const noexcept;

private:
    static unsigned int newUID;
    int total_alleles_;
};


#endif //TRANSMISSION_NETWORKS_APP_LOCUS_H
