//
// Created by Maxwell Murphy on 10/16/19.
//

#ifndef MALARIA_NET_EDGE_H
#define MALARIA_NET_EDGE_H

#include "Node.h"

struct Edge {
    int num_generations;
    Node& parent;
    Node& child;
};

#endif //MALARIA_NET_EDGE_H
