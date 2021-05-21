//
// Created by Maxwell Murphy on 2/25/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_FORWARDING_UTILS_H
#define TRANSMISSION_NETWORKS_APP_FORWARDING_UTILS_H


// Technique to deal with perfect forwarding problems
// https://akrzemi1.wordpress.com/2013/10/10/too-perfect-forwarding/

#define ENABLE_IF(...)  typename std::enable_if<__VA_ARGS__>::type* = nullptr

namespace transmission_nets::core::utils {

    template <typename T, typename U>
    constexpr inline bool NonSelf()
    {
        using DecayedT = typename std::decay<T>::type;
        return !std::is_same<DecayedT, U>::value
               && !std::is_base_of<U, DecayedT>::value;
    }

}


#endif //TRANSMISSION_NETWORKS_APP_FORWARDING_UTILS_H
