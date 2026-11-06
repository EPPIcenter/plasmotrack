//
// A prior term given as a function of its inputs: the log density `f()` reads, recomputed when any input changes.
// For priors with no class of their own (the transmission leak's, the hierarchical observation rates').
//

#ifndef TRANSMISSION_NETWORKS_APP_FUNCTIONLOGPDF_H
#define TRANSMISSION_NETWORKS_APP_FUNCTIONLOGPDF_H

#include "core/abstract/observables/Dependency.h"
#include "core/computation/PartialLikelihood.h"

#include <functional>
#include <string>
#include <utility>

namespace transmission_nets::core::distributions {

    class FunctionLogPDF : public computation::PartialLikelihood {
    public:
        template<typename... Inputs>
        FunctionLogPDF(std::string id, std::function<computation::Likelihood()> f, Inputs&&... inputs) : id_(std::move(id)), f_(std::move(f)) {
            abstract::dependsOn(this, std::forward<Inputs>(inputs)...);
            this->setDirty();
        }

        computation::Likelihood value() override {
            if (isDirty()) {
                value_ = f_();
                this->setClean();
            }
            return value_;
        }

        std::string identifier() override { return id_; }

    private:
        std::string id_;
        std::function<computation::Likelihood()> f_;
    };

}// namespace transmission_nets::core::distributions

#endif//TRANSMISSION_NETWORKS_APP_FUNCTIONLOGPDF_H
