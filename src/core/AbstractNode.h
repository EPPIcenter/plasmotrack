//
// Created by Maxwell Murphy on 12/30/19.
//

#ifndef TRANSMISSION_NETWORKS_APP_ABSTRACTNODE_H
#define TRANSMISSION_NETWORKS_APP_ABSTRACTNODE_H

#include <boost/container/flat_set.hpp>

class AbstractNode {
public:
    virtual bool isDirty() = 0;
    virtual void setDirty() = 0;
    virtual void setDirty(AbstractNode* ptr) = 0;
    virtual void setClean() = 0;
    virtual void saveState() = 0;
    virtual void restoreState() = 0;
    virtual void acceptState() = 0;
    virtual void addDependency(AbstractNode* ptr) = 0;
};


#endif //TRANSMISSION_NETWORKS_APP_ABSTRACTNODE_H
