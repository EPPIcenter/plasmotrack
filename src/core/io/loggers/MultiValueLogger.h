//
// Created by Maxwell Murphy on 5/18/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_MULTIVALUELOGGER_H
#define TRANSMISSION_NETWORKS_APP_MULTIVALUELOGGER_H

#include "AbstractLogger.h"
#include "core/io/serialize.h"
#include <string>

namespace transmission_nets::core::io {

    template<typename T>
    class MultiValueLogger : public AbstractLogger {

    public:
        template<typename Output>
        explicit MultiValueLogger(std::unique_ptr<Output> output);

        std::string prepareValue() noexcept override;

        void addTarget(const std::string& id, std::shared_ptr<T> target);


    private:
        std::vector<std::pair<std::string, std::shared_ptr<T>>> targets_{};
        int iter_ = 1;
    };

    template<typename T>
    template<typename Output>
    MultiValueLogger<T>::MultiValueLogger(std::unique_ptr<Output> output) : AbstractLogger(std::move(output)) {}

    template<typename T>
    std::string MultiValueLogger<T>::prepareValue() noexcept {
        std::stringstream ss;
        for (const auto& [id, target] : targets_) {
            ss << id << "," << serialize(target->value()) << "," << iter_ << "\n";
        }
        ++iter_;
        return ss.str();
    }

    template<typename T>
    void MultiValueLogger<T>::addTarget(const std::string& id, std::shared_ptr<T> target) {
        targets_.push_back({id, target});
    }


}// namespace transmission_nets::core::io

#endif//TRANSMISSION_NETWORKS_APP_MULTIVALUELOGGER_H
