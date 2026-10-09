//
// Groups infections by latent genome, for samplers whose proposals depend on other infections only through them.
//

#ifndef TRANSMISSION_NETWORKS_APP_GENOMEPARTITION_H
#define TRANSMISSION_NETWORKS_APP_GENOMEPARTITION_H

#include <algorithm>
#include <cstdint>
#include <vector>

namespace transmission_nets::core::samplers::genetics {

    /*
     * Partitions infections into classes with identical latent genomes. A proposal that depends on other infections
     * only through their genomes needs evaluating once per class; in a clonal cluster most candidate parents share one.
     *
     * Samplers evaluate their proposal before and after changing only their own infection, so given the same
     * infections with unchanged genomes, the partition is kept rather than rebuilt.
     */
    template<typename InfectionEventImpl>
    class GenomePartition {
    public:
        // Partitions a range of pointers to infections. Classes are numbered in order of first appearance.
        template<typename Range>
        void update(const Range& infections);

        [[nodiscard]] size_t numClasses() const noexcept { return members_.size(); }
        // The class of the i-th infection given to update().
        [[nodiscard]] std::uint32_t classOf(const size_t i) const noexcept { return classOf_[i]; }
        [[nodiscard]] std::uint32_t size(const std::uint32_t c) const noexcept { return sizes_[c]; }
        // A member of the class, whose genome is the class's.
        [[nodiscard]] const InfectionEventImpl& member(const std::uint32_t c) const noexcept { return *members_[c]; }

        // Whether two infections have the same latent genotype at every locus.
        static bool sameGenome(const InfectionEventImpl& a, const InfectionEventImpl& b);

    private:
        std::vector<const InfectionEventImpl*> infections_{};
        std::vector<std::uint64_t> versions_{};
        std::vector<std::uint32_t> classOf_{};
        std::vector<std::uint32_t> sizes_{};
        std::vector<const InfectionEventImpl*> members_{};
    };

    template<typename InfectionEventImpl>
    template<typename Range>
    void GenomePartition<InfectionEventImpl>::update(const Range& infections) {
        size_t n       = 0;
        bool unchanged = true;
        for (const auto& infection : infections) {
            unchanged = unchanged and n < infections_.size() and infections_[n] == &*infection and versions_[n] == infection->latentGenomeVersion();
            ++n;
        }
        if (unchanged and n == infections_.size()) {
            return;
        }

        infections_.clear();
        versions_.clear();
        classOf_.clear();
        sizes_.clear();
        members_.clear();
        for (const auto& infection : infections) {
            const InfectionEventImpl* added = &*infection;
            infections_.push_back(added);
            versions_.push_back(added->latentGenomeVersion());
            std::uint32_t c = 0;
            while (c < members_.size() and not sameGenome(*members_[c], *added)) {
                ++c;
            }
            if (c == members_.size()) {
                members_.push_back(added);
                sizes_.push_back(0);
            }
            ++sizes_[c];
            classOf_.push_back(c);
        }
    }

    template<typename InfectionEventImpl>
    bool GenomePartition<InfectionEventImpl>::sameGenome(const InfectionEventImpl& a, const InfectionEventImpl& b) {
        // Genotypes are kept by locus in the same order in every infection, so a single pass compares them.
        const auto& first  = a.latentGenotype();
        const auto& second = b.latentGenotype();
        return first.size() == second.size() and std::equal(first.begin(), first.end(), second.begin(), [](const auto& x, const auto& y) {
                   return x.first == y.first and x.second->value() == y.second->value();
               });
    }

}// namespace transmission_nets::core::samplers::genetics

#endif//TRANSMISSION_NETWORKS_APP_GENOMEPARTITION_H
