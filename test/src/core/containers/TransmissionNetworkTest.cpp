//
// Created by Maxwell Murphy on 6/10/20.
//

#include "gtest/gtest.h"


#include "core/containers/Infection.h"
#include "core/datatypes/Alleles.h"
#include "core/parameters/TransmissionNetwork.h"

using namespace transmission_nets::core::datatypes;
using namespace transmission_nets::core::containers;
using namespace transmission_nets::core::parameters;

TEST(TransmissionNetworkTest, CoreTest) {
    using NodeImpl = Infection<AllelesBitSet<16>>;
    TransmissionNetwork<NodeImpl> network;

    bool inf1ParentSetChanged = false;
    bool inf2ParentSetChanged = false;
    bool inf3ParentSetChanged = false;
    bool inf4ParentSetChanged = false;

    auto inf1 = new NodeImpl("node1", 1);
    auto inf2 = new NodeImpl("node2", 1);
    auto inf3 = new NodeImpl("node3", 1);
    auto inf4 = new NodeImpl("node4", 1);


    network.addNode(inf1);
    network.addNode(inf2);
    network.addNode(inf3);
    network.addNode(inf4);

    ASSERT_EQ(network.nodes().size(), 4);

    network.parentSet(inf1)->add_post_change_listener([&]() { inf1ParentSetChanged = true; });
    network.parentSet(inf2)->add_post_change_listener([&]() { inf2ParentSetChanged = true; });
    network.parentSet(inf3)->add_post_change_listener([&]() { inf3ParentSetChanged = true; });
    network.parentSet(inf4)->add_post_change_listener([&]() { inf4ParentSetChanged = true; });

    ASSERT_FALSE(network.createsCycle(inf4, inf1));
    ASSERT_FALSE(network.createsCycle(inf3, inf1));
    ASSERT_FALSE(network.createsCycle(inf2, inf1));

    ASSERT_FALSE(inf1ParentSetChanged);
    ASSERT_FALSE(inf2ParentSetChanged);
    ASSERT_FALSE(inf3ParentSetChanged);
    ASSERT_FALSE(inf4ParentSetChanged);

    network.parentSet(inf2)->saveState("state1");
    network.addEdge(inf1, inf2);
    network.parentSet(inf2)->acceptState();

    network.parentSet(inf3)->saveState("state1");
    network.addEdge(inf2, inf3);
    network.parentSet(inf3)->acceptState();

    network.parentSet(inf4)->saveState("state1");
    network.addEdge(inf3, inf4);
    network.parentSet(inf4)->acceptState();


    ASSERT_FALSE(inf1ParentSetChanged);
    ASSERT_TRUE(inf2ParentSetChanged);
    ASSERT_TRUE(inf3ParentSetChanged);
    ASSERT_TRUE(inf4ParentSetChanged);

    ASSERT_TRUE(network.createsCycle(inf4, inf1));
    ASSERT_TRUE(network.createsCycle(inf3, inf1));
    ASSERT_TRUE(network.createsCycle(inf2, inf1));
    ASSERT_TRUE(network.createsCycle(inf1, inf1));

}