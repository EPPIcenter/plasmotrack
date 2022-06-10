//
// Created by Maxwell Murphy on 12/9/19.
//

#ifndef TRANSMISSION_NETWORKS_APP_CRTP_H
#define TRANSMISSION_NETWORKS_APP_CRTP_H

// Helper for the CRTP design pattern
// https://www.fluentcpp.com/2017/05/19/crtp-helper/

namespace transmission_nets::core::abstract {

    template<typename T, template<typename...> typename crtpType, typename... Us>
    struct crtp {
        T& underlying() { return static_cast<T&>(*this); }
        T const& underlying() const { return static_cast<T const&>(*this); }

    private:
        crtp() = default;
        friend crtpType<T, Us...>;
    };

}// namespace transmission_nets::core::abstract


#endif//TRANSMISSION_NETWORKS_APP_CRTP_H
