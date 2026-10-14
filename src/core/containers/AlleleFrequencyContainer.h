//
// Created by Maxwell Murphy on 3/6/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCYCONTAINER_H
#define TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCYCONTAINER_H

#include "core/abstract/observables/CheckpointablePassthrough.h"
#include "core/abstract/observables/Dependency.h"
#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/UncacheablePassthrough.h"

#include "core/containers/Locus.h"
#include "core/containers/UidLess.h"

#include "core/parameters/Parameter.h"

#include <boost/container/flat_map.hpp>

#include <memory>
#include <unordered_map>


namespace transmission_nets::core::containers {

    template<typename AlleleFrequencyImpl, typename LocusImpl = core::containers::Locus>
    class AlleleFrequencyContainer : public abstract::Observable<AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>>,
                                     public abstract::UncacheablePassthrough<AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>>,
                                     public abstract::CheckpointablePassthrough<AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>> {

        using AlleleFrequencyMap = boost::container::flat_map<std::shared_ptr<LocusImpl>, std::shared_ptr<parameters::Parameter<AlleleFrequencyImpl>>, UidLess>;

    public:
        using LocusAlleleFrequencyAssignment = std::pair<std::shared_ptr<LocusImpl>, std::shared_ptr<AlleleFrequencyImpl>>;

        std::vector<std::shared_ptr<LocusImpl>> loci{};

        AlleleFrequencyContainer() = default;


        std::shared_ptr<parameters::Parameter<AlleleFrequencyImpl>> alleleFrequencies(std::shared_ptr<LocusImpl> locus) {
            return alleleFrequencies_.at(locus);
        };

        void addLocus(std::shared_ptr<LocusImpl> locus);

        double totalLoci() {
            return alleleFrequencies_.size();
        }

        friend std::ostream& operator<<(std::ostream& os, const AlleleFrequencyContainer& container) {
            for (const auto& [locus, alleleFreqs] : container.alleleFrequencies_) {
                os << locus->label << ": " << alleleFreqs.value() << "\n";
            }
            return os;
        };


    private:
        AlleleFrequencyMap alleleFrequencies_{};
    };

    template<typename AlleleFrequencyImpl, typename LocusImpl>
    void AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>::addLocus(std::shared_ptr<LocusImpl> locus) {
        auto af = std::make_shared<parameters::Parameter<AlleleFrequencyImpl>>(locus->totalAlleles());
        loci.push_back(locus);
        alleleFrequencies_.insert_or_assign(locus, af);
        abstract::passThrough(this, af);
    }


}// namespace transmission_nets::core::containers


#endif//TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCYCONTAINER_H
