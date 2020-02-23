//
// Created by Maxwell Murphy on 1/23/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_LOGLIKELIHOOD_H
#define TRANSMISSION_NETWORKS_APP_LOGLIKELIHOOD_H

#include <functional>
#include <optional>
#include <utility>
#include <boost/container/flat_map.hpp>


class LogLikelihood {
    using Callback = std::function<void(const LogLikelihood)>;
    template<typename Key, typename Value> using Map_ = boost::container::flat_map<Key, Value>;

public:
    enum callbackId: uint_fast64_t {};

    LogLikelihood(std::string id);

    LogLikelihood(float value, std::string id);

    float value();

private:
    float value_;
    std::string id_;
    std::optional<float> saved_state_{};

    Map_<uint_fast64_t, Callback> pre_change_callbacks_{};
    Map_<uint_fast64_t, Callback> post_change_callbacks_{};
    Map_<uint_fast64_t, Callback> save_state_callbacks_{};
    Map_<uint_fast64_t, Callback> accept_state_callbacks_{};
    Map_<uint_fast64_t, Callback> restore_state_callbacks_{};
};


#endif //TRANSMISSION_NETWORKS_APP_LOGLIKELIHOOD_H
