//
// Created by Maxwell Murphy on 12/4/19.
//

#ifndef TRANSMISSION_NETWORKS_CORE_PARAMETER_H
#define TRANSMISSION_NETWORKS_CORE_PARAMETER_H

#include <string>
#include "../Node.h"

// A node that:
// 1: has no parent dependencies
// 2: may have its value explicitly set
// 3: may have its value updated
// Serves as an input to the computational graph


template<typename T>
class Parameter : public Node<T> {
protected:
    std::string id_;

public:
    Parameter(const std::string& id, T value);
    Parameter(const std::string& id);

    T value() override;

    T peek() override;

    void setValue(T val);

    std::string id();

};

template<typename T>
void Parameter<T>::setValue(T value) {
    this->setDirty();
    this->value_ = value;
    this->setClean();
}

template<typename T>
Parameter<T>::Parameter(const std::string& id, T value) {
    this->id_ = id;
    this->value_ = value;
    this->setClean();
}

template<typename T>
Parameter<T>::Parameter(const std::string& id) {
    this->id_ = id;
}

template<typename T>
T Parameter<T>::value() {
    return this->value_;
}

template<typename T>
T Parameter<T>::peek() {
    return this->value_;
}

template<typename T>
std::string Parameter<T>::id() {
    return this->id_;
}


#endif //TRANSMISSION_NETWORKS_CORE_PARAMETER_H
