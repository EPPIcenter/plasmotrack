//
// Created by Maxwell Murphy on 2/6/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_GRAPH_H
#define TRANSMISSION_NETWORKS_APP_GRAPH_H

#include <boost/container/flat_map.hpp>

#include "core/containers/ParentSet.h"


template<typename NodeValueImpl>
class TransmissionNetwork :  public Observable<TransmissionNetwork<NodeValueImpl>>,
                             public UncacheablePassthrough<TransmissionNetwork<NodeValueImpl>>,
                             public CheckpointablePassthrough<TransmissionNetwork<NodeValueImpl>> {


public:

private:
    boost::container::flat_map<NodeValueImpl*, Parameter<ParentSet<NodeValueImpl>>> parentSets_{};

};


#endif //TRANSMISSION_NETWORKS_APP_GRAPH_H
