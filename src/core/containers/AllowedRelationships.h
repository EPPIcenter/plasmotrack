//
// Created by Maxwell Murphy on 1/10/24.
//

#ifndef ALLOWEDRELATIONSHIPS_H
#define ALLOWEDRELATIONSHIPS_H

#include <map>
#include <memory>
#include <vector>


namespace transmission_nets::core::containers {
    template<typename InfectionEvent>
    struct AllowedRelationships {
        std::map<std::shared_ptr<InfectionEvent>, std::vector<std::shared_ptr<InfectionEvent>>> allowedParents;
        std::map<std::shared_ptr<InfectionEvent>, std::vector<std::shared_ptr<InfectionEvent>>> allowedChildren;
    };
}


#endif //ALLOWEDRELATIONSHIPS_H
