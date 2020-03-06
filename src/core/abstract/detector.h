//
// Created by Maxwell Murphy on 2/28/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_DETECTOR_H
#define TRANSMISSION_NETWORKS_APP_DETECTOR_H

template <template <class> class Op, class T, class = void>
struct is_detected : std::false_type {};
template <template <class> class Op, class T>
struct is_detected<Op, T, std::void_t<Op<T>>> : std::true_type {};

template <template <class> class Op, class T>
inline constexpr bool is_detected_v = is_detected<Op, T>::value;

#endif //TRANSMISSION_NETWORKS_APP_DETECTOR_H
