//
// Created by Maxwell Murphy on 2/21/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESSV3_H
#define TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESSV3_H

#include "core/computation/OrderDerivedParentSet.h"
#include "core/computation/PartialLikelihood.h"
#include "core/containers/Infection.h"
#include "core/io/serialize.h"
#include "core/utils/generators/CombinationIndicesGenerator.h"
#include "core/utils/numerics.h"


#include <boost/container/flat_set.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/range/adaptors.hpp>
#include <fmt/core.h>
//#include <fmt/ranges.h>

namespace transmission_nets::model::transmission_process {

    using Likelihood = core::computation::Likelihood;

    template<typename InfectionEventImpl>
    struct ParentSetDist {
        Likelihood sourceLlik = 0;
        std::vector<std::pair<Likelihood, core::containers::ParentSet<InfectionEventImpl>>> parentSetLliks{};
        Likelihood totalLlik = 0;
    };

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    class OrderBasedTransmissionProcessV3 : public core::computation::PartialLikelihood {

        /*
         * The order based transmission process considers the set of all possible parent sets under the given ordering.
         * Assumes each infection has a potential latent unobserved parent that is drawn from the background distribution.
         */

        using ListenerIdMap = boost::container::flat_map<std::string, core::abstract::ListenerId_t>;

    public:
        OrderBasedTransmissionProcessV3(std::shared_ptr<NodeTransmissionProcessImpl> ntp,
                                        std::shared_ptr<SourceTransmissionProcessImpl> stp,
                                        std::shared_ptr<InfectionEventImpl> child,
                                        std::shared_ptr<ParentSetImpl> parent_set,
                                        std::shared_ptr<InfectionEventImpl> latent_parent);

        Likelihood value() override;

        // Allows querying the individual parent set likelihoods
        ParentSetDist<InfectionEventImpl> calcParentSetDist();

        std::string identifier() override;

        Likelihood peek() noexcept override;

        // Input parameters
        std::shared_ptr<NodeTransmissionProcessImpl> ntp_;
        std::shared_ptr<SourceTransmissionProcessImpl> stp_;
        std::shared_ptr<InfectionEventImpl> child_;
        std::shared_ptr<ParentSetImpl> parentSet_;
        std::shared_ptr<InfectionEventImpl> latentParent_;

    private:
        void postSaveState(const std::string& savedStateId);

        void postAcceptState();

        void postRestoreState(const std::string& savedStateId);

        Likelihood getLikelihood(const core::containers::ParentSet<InfectionEventImpl> &ps);
        bool likelihoodCalculated(const core::containers::ParentSet<InfectionEventImpl> &ps);
        void setLikelihood(const core::containers::ParentSet<InfectionEventImpl> &ps, Likelihood llik);
        void clearParentLikelihood(std::shared_ptr<InfectionEventImpl> parent);
        void clearLikelihood();

        // Container to track the calculated parent set likelihoods over which we sum
        // to get the total likelihood
        // todo: figure out how to reduce overhead of tracking. currently violates state tracking...
        using LikelihoodTracker = boost::container::flat_map<boost::container::flat_set<std::string>, Likelihood>;

        std::array<LikelihoodTracker, 25> parentSetLliks_{};
        size_t parentSetLLiksIndex_ = 0;

//        LikelihoodTracker *parentSetLliks_ = new LikelihoodTracker {};
//        std::deque<LikelihoodTracker*> parentSetLliksCache_{};

    };

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::OrderBasedTransmissionProcessV3(
            std::shared_ptr<NodeTransmissionProcessImpl> ntp,
            std::shared_ptr<SourceTransmissionProcessImpl> stp,
            std::shared_ptr<InfectionEventImpl> child,
            std::shared_ptr<ParentSetImpl> parent_set,
            std::shared_ptr<InfectionEventImpl> latent_parent) : ntp_(std::move(ntp)),
                                                                 stp_(std::move(stp)),
                                                                 child_(std::move(child)),
                                                                 parentSet_(std::move(parent_set)),
                                                                 latentParent_(std::move(latent_parent)) {

        ntp_->add_set_dirty_listener([=, this]() {
            clearLikelihood();
            setDirty();
        });
        ntp_->registerCacheableCheckpointTarget(this);

        stp_->add_set_dirty_listener([=, this]() {
            clearParentLikelihood(latentParent_);
            setDirty();
        });
        stp_->registerCacheableCheckpointTarget(this);

        child_->add_post_change_listener([=, this]() {
            clearLikelihood();
            setDirty();
        });
        child_->registerCacheableCheckpointTarget(this);

        latentParent_->add_post_change_listener([=, this]() {
            clearParentLikelihood(latentParent_);
            setDirty();
        });
        latentParent_->registerCacheableCheckpointTarget(this);

        parentSet_->add_element_changed_listener([=, this](const auto &parent) {
            clearParentLikelihood(parent);
            setDirty();
        });

        parentSet_->add_element_added_listener([=, this]([[maybe_unused]] const auto &parent) {
            setDirty();
        });

        parentSet_->add_element_removed_listener([=, this](const auto &parent) {
            clearParentLikelihood(parent);
            setDirty();
        });

        parentSet_->registerCacheableCheckpointTarget(this);

        this->addPostSaveHook([=, this](const auto savedStateID) { this->postSaveState(savedStateID); });
        this->addPostRestoreHook([=, this](const auto savedStateID) { this->postRestoreState(savedStateID); });
        this->addPostAcceptHook([=, this]() { this->postAcceptState(); });

        this->value_ = -std::numeric_limits<Likelihood>::infinity();
        this->setDirty();
        this->value();

#ifndef NDEBUG
        if (this->value_ <= -std::numeric_limits<Likelihood>::infinity()) {
            fmt::print("Child: {}\n", this->child_->id());
            fmt::print("Parents: ");
            for (const auto& parent : this->parentSet_->value()) {
                fmt::print("{}, ", parent->id());
            }
            fmt::print("\n");
            fmt::print("LLiks: \n");
            for (const auto& [ps, llik]: parentSetLliks_) {
                fmt::print("\t{}, {}\n", core::io::serialize(ps), llik);
            }
            fmt::print("\n");
            exit(1);
        }
#endif
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    Likelihood OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::getLikelihood(const core::containers::ParentSet<InfectionEventImpl> &ps) {
        boost::container::flat_set<std::string> key;
        for (const auto &parent : ps) {
            key.insert(parent->id());
        }
#ifndef NDEBUG
        if (!parentSetLliks_.contains(key)) {
            fmt::print(stderr, "Likelihood not found for parent set: ");
            for (const auto &parent : ps) {
                fmt::print(stderr, "{} ", parent->id());
            }
            fmt::print(stderr, "\n");
            fmt::print(stderr, "Valid parent sets: \n");
            for (const auto &parentSet : parentSetLliks_) {
                fmt::print(stderr, "\tParent set: ");
                for (const auto &parent : parentSet.first) {
                    fmt::print(stderr, "{} ", parent);
                }
                fmt::print(stderr, "\n");
            }

            exit(1);
        }
#endif
        return parentSetLliks_[parentSetLLiksIndex_].at(key);
    }

    //likelihood calculated
    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    bool OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::likelihoodCalculated(const core::containers::ParentSet<InfectionEventImpl> &ps) {
        boost::container::flat_set<std::string> key;
        for (const auto &parent : ps) {
            key.insert(parent->id());
        }
        return parentSetLliks_[parentSetLLiksIndex_].contains(key);
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::setLikelihood(const core::containers::ParentSet<InfectionEventImpl> &ps, Likelihood llik) {
        boost::container::flat_set<std::string> key;
        for (const auto &parent : ps) {
            key.insert(parent->id());
        }
        parentSetLliks_[parentSetLLiksIndex_][key] = llik;
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::clearParentLikelihood(const std::shared_ptr<InfectionEventImpl> parent) {
        for (auto it = parentSetLliks_[parentSetLLiksIndex_].begin(); it != parentSetLliks_[parentSetLLiksIndex_].end();) {
            if (it->first.find(parent->id()) != it->first.end()) {
                it = parentSetLliks_[parentSetLLiksIndex_].erase(it);
            } else {
                ++it;
            }
        }
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::clearLikelihood() {
        parentSetLliks_[parentSetLLiksIndex_].clear();
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    std::string OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::identifier() {
        auto out = fmt::format("OrderBasedTransmissionProcessV3<{}>", child_->id());
        return out;
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    Likelihood
    OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::value() {
        if (this->isDirty()) {

            core::utils::generators::CombinationIndicesGenerator comboGen;
            std::vector<Likelihood> lliks{};
            Likelihood ps_llik;
            Likelihood maxLlik = -std::numeric_limits<Likelihood>::infinity();

            const auto ps = parentSet_->value();
            const int totalNodes = ps.size();

            // Calculate the single latent parent case
            core::containers::ParentSet<InfectionEventImpl> tmpPs_{latentParent_};

            if (likelihoodCalculated(tmpPs_)) {
                ps_llik = getLikelihood(tmpPs_);
            } else {
                ps_llik = ntp_->calculateLogLikelihood(child_, latentParent_, stp_);
                setLikelihood(tmpPs_, ps_llik);
            }
            lliks.push_back(ps_llik);
            maxLlik = std::max(maxLlik, lliks.back());

            // Iterate over all possible parent sets
            for (int i = 1; i <= ParentSetMaxCardinality and i <= totalNodes; i++) {
                comboGen.reset(totalNodes, i);
                while (!comboGen.completed) {
                    // Generate the parent set
                    tmpPs_.clear();
                    for (const auto& idx : comboGen.curr) {
                        tmpPs_.insert(ps.begin()[idx]);
                    }

                    // Calculate the likelihood without latent parent
                    if (likelihoodCalculated(tmpPs_)) {
                        ps_llik = getLikelihood(tmpPs_);
                    } else {
                        ps_llik = ntp_->calculateLogLikelihood(child_, tmpPs_);
                        setLikelihood(tmpPs_, ps_llik);
                    }
                    lliks.push_back(ps_llik);
                    maxLlik = std::max(maxLlik, lliks.back());

                    // Calculate with latent parent

                    tmpPs_.insert(latentParent_);
                    if (likelihoodCalculated(tmpPs_)) {
                        ps_llik = getLikelihood(tmpPs_);
                    } else {
                        tmpPs_.erase(latentParent_);
                        ps_llik = ntp_->calculateLogLikelihood(child_, latentParent_, tmpPs_, stp_);
                        tmpPs_.insert(latentParent_);
                        setLikelihood(tmpPs_, ps_llik);
                    }
                    lliks.push_back(ps_llik);
                    maxLlik = std::max(maxLlik, lliks.back());

                    comboGen.next();
                }
            }

            this->value_ = core::utils::logSumExpKnownMax(lliks.begin(), lliks.end(), maxLlik);

            assert(this->value_ < std::numeric_limits<Likelihood>::infinity());

            this->setClean();
        }
        return this->value_;
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::postSaveState([[maybe_unused]] const std::string& savedStateId) {
        parentSetLliks_[parentSetLLiksIndex_ + 1] = parentSetLliks_[parentSetLLiksIndex_];
        parentSetLLiksIndex_++;
//        parentSetLliksCache_.push_back(parentSetLliks_);
//        parentSetLliks_ = new LikelihoodTracker {*parentSetLliksCache_.back()};
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::postRestoreState([[maybe_unused]] const std::string& savedStateId) {
        parentSetLLiksIndex_--;
//        delete parentSetLliks_;
//        parentSetLliks_ = parentSetLliksCache_.back();
//        parentSetLliksCache_.pop_back();
    }


    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    void OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::postAcceptState() {
        parentSetLliks_[0] = parentSetLliks_[parentSetLLiksIndex_];
        parentSetLLiksIndex_ = 0;
//        for (auto el : parentSetLliksCache_) {
//            delete el;
//        }
//        parentSetLliksCache_.clear();
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    Likelihood OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::peek() noexcept {

#ifndef NDEBUG
        if (this->value_ <= -std::numeric_limits<Likelihood>::infinity()) {
            std::cerr << "Saved States: " << this->saved_states_stack_.size() << std::endl;
            for ([[maybe_unused]] const auto& savedState : this->saved_states_stack_) {
                std::cerr << "State: " << savedState.saved_state << " " << savedState.saved_state_id << std::endl;
            }
        }
#endif
        return Computation::peek();
    }

    template<int ParentSetMaxCardinality, typename NodeTransmissionProcessImpl, typename SourceTransmissionProcessImpl, typename InfectionEventImpl, typename ParentSetImpl>
    ParentSetDist<InfectionEventImpl> OrderBasedTransmissionProcessV3<ParentSetMaxCardinality, NodeTransmissionProcessImpl, SourceTransmissionProcessImpl, InfectionEventImpl, ParentSetImpl>::calcParentSetDist() {
        ParentSetDist<InfectionEventImpl> dist{};
        dist.totalLlik = this->value();
        const auto tmpPs = parentSet_->value();
        const int totalNodes = tmpPs.size();
        core::containers::ParentSet<InfectionEventImpl> ps{};
        ps.insert(latentParent_);
        dist.parentSetLliks.push_back(std::make_pair(getLikelihood(ps), ps));

        core::utils::generators::CombinationIndicesGenerator comboGen;
        for (int i = 1; i <= ParentSetMaxCardinality; i++) {
            comboGen.reset(totalNodes, i);
            while (!comboGen.completed) {
                ps.clear();
                for (const auto& idx : comboGen.curr) {
                    ps.insert(tmpPs.begin()[idx]);
                }

//              calculate parentset without latent parent
                dist.parentSetLliks.push_back(std::make_pair(getLikelihood(ps), ps));

//              calculate parentset with latent parent
                ps.insert(latentParent_);
                dist.parentSetLliks.push_back(std::make_pair(getLikelihood(ps), ps));
                comboGen.next();
            }
        }

        return dist;
    }

}// namespace transmission_nets::model::transmission_process


#endif//TRANSMISSION_NETWORKS_APP_ORDERBASEDTRANSMISSIONPROCESSV3_H
