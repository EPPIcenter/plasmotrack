//
// Created by Maxwell Murphy on 11/25/19.
//

#ifndef TRANSMISSION_NETWORKS_APP_MODELPARAMETER_H
#define TRANSMISSION_NETWORKS_APP_MODELPARAMETER_H

template <class T>
class ModelParameter {
private:
    std::vector<std::shared_ptr<PartialLikelihood>> referring_pliks;
public:
    T previous_value;
    T value;
    virtual void update() = 0;
    void register_plik(std::shared_ptr<PartialLikelihood> plik);
};


#endif //TRANSMISSION_NETWORKS_APP_MODELPARAMETER_H
