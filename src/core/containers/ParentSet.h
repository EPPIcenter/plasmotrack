//
// Created by Maxwell Murphy on 2/19/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_PARENTSET_H
#define TRANSMISSION_NETWORKS_APP_PARENTSET_H

#include <boost/container/flat_set.hpp>

#include <memory>


namespace transmission_nets::core::containers {

    template<typename ElementType>
    using ParentSet = boost::container::flat_set<std::shared_ptr<ElementType>>;

}

#endif //TRANSMISSION_NETWORKS_APP_PARENTSET_H
