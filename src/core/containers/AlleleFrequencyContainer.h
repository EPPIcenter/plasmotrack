//
// Created by Maxwell Murphy on 3/6/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCYCONTAINER_H
#define TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCYCONTAINER_H

#include "core/abstract/observables/Observable.h"
#include "core/abstract/observables/UncacheablePassthrough.h"
#include "core/abstract/observables/CheckpointablePassthrough.h"

#include "core/containers/Locus.h"

#include "core/parameters/Parameter.h"

#include <memory>
#include <unordered_map>


namespace transmission_nets::core::containers {

    template<typename AlleleFrequencyImpl, typename LocusImpl = core::containers::Locus>
    class AlleleFrequencyContainer : public abstract::Observable<AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>>,
                                     public abstract::UncacheablePassthrough<AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>>,
                                     public abstract::CheckpointablePassthrough<AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>> {

        using AlleleFrequencyMap = std::unordered_map<std::shared_ptr<LocusImpl>, std::shared_ptr<parameters::Parameter<AlleleFrequencyImpl>>>;

    public:
        using LocusAlleleFrequencyAssignment = std::pair<std::shared_ptr<LocusImpl>, std::shared_ptr<AlleleFrequencyImpl>>;

        std::vector<std::shared_ptr<LocusImpl>> loci{};

        AlleleFrequencyContainer() = default;

//        template<typename AlleleDataIter>
//        explicit AlleleFrequencyContainer(AlleleDataIter alleleFreqs);

        std::shared_ptr<parameters::Parameter<AlleleFrequencyImpl>> alleleFrequencies(std::shared_ptr<LocusImpl> locus) {
            return alleleFrequencies_.at(locus);
        };

        void addLocus(std::shared_ptr<LocusImpl> locus);

        double totalLoci() {
            return alleleFrequencies_.size();
        }

        friend std::ostream &operator<<(std::ostream &os, const AlleleFrequencyContainer &container) {
            for (const auto& [locus, alleleFreqs] : container.alleleFrequencies_ ) {
                os << locus->label << ": " << alleleFreqs.value() << "\n";
            }
            return os;
        };


    private:
        AlleleFrequencyMap alleleFrequencies_{};

    };

//    template<typename AlleleFrequencyImpl, typename LocusImpl>
//    template<typename AlleleDataIter>
//    AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>::AlleleFrequencyContainer(AlleleDataIter alleleFreqs) {
//        for (auto &[locus, alleleFrequencies] : alleleFreqs) {
//            assert(locus->totalAlleles() == alleleFrequencies.totalElements());
//            auto af = std::make_shared<parameters::Parameter<AlleleFrequencyImpl>>(*alleleFrequencies);
//            loci.push_back(locus);
//            alleleFrequencies_.insert_or_assign(locus, af);
//            // pass through notifications
//            alleleFrequencies_.at(locus)->add_pre_change_listener([=, this]() { this->notify_pre_change(); });
//            alleleFrequencies_.at(locus)->add_post_change_listener([=, this]() { this->notify_post_change(); });
//            alleleFrequencies_.at(locus)->add_save_state_listener([=, this](const std::string saveStateId) { this->notify_save_state(saveStateId); });
//            alleleFrequencies_.at(locus)->add_accept_state_listener([=, this]() { this->notify_accept_state(); });
//            alleleFrequencies_.at(locus)->add_restore_state_listener([=, this](const std::string saveStateId) { this->notify_restore_state(saveStateId); });
//        }
//    }

    template<typename AlleleFrequencyImpl, typename LocusImpl>
    void AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>::addLocus(std::shared_ptr<LocusImpl> locus) {
        auto af = std::make_shared<parameters::Parameter<AlleleFrequencyImpl>>(locus->totalAlleles());
        loci.push_back(locus);
        alleleFrequencies_.insert_or_assign(locus, af);
        alleleFrequencies_.at(locus)->add_pre_change_listener([=, this]() { this->notify_pre_change(); });
        alleleFrequencies_.at(locus)->add_post_change_listener([=, this]() { this->notify_post_change(); });
        alleleFrequencies_.at(locus)->add_save_state_listener([=, this](const std::string saveStateId) { this->notify_save_state(saveStateId); });
        alleleFrequencies_.at(locus)->add_accept_state_listener([=, this]() { this->notify_accept_state(); });
        alleleFrequencies_.at(locus)->add_restore_state_listener([=, this](const std::string saveStateId) { this->notify_restore_state(saveStateId); });
    }


}


#endif //TRANSMISSION_NETWORKS_APP_ALLELEFREQUENCYCONTAINER_H
