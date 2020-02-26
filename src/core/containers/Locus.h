//
// Created by Maxwell Murphy on 2/23/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_LOCUS_H
#define TRANSMISSION_NETWORKS_APP_LOCUS_H

class Locus {
public:
    const unsigned int uid;
    const std::string label;

    explicit Locus(std::string label);

    bool operator<(const Locus &rhs) const;

    bool operator>(const Locus &rhs) const;

    bool operator<=(const Locus &rhs) const;

    bool operator>=(const Locus &rhs) const;

private:
    static unsigned int newUID;
};


#endif //TRANSMISSION_NETWORKS_APP_LOCUS_H
