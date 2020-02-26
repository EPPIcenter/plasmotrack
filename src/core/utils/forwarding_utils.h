//
// Created by Maxwell Murphy on 2/25/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_FORWARDING_UTILS_H
#define TRANSMISSION_NETWORKS_APP_FORWARDING_UTILS_H

#define ENABLE_IF(...) \
  typename std::enable_if<__VA_ARGS__>::type* = nullptr

template <typename T, typename U>
constexpr inline bool NonSelf()
{
    using DecayedT = typename std::decay<T>::type;
    return !std::is_same<DecayedT, U>::value
           && !std::is_base_of<U, DecayedT>::value;
}


#endif //TRANSMISSION_NETWORKS_APP_FORWARDING_UTILS_H
