//
// Version information for Plasmotrack
//

#ifndef TRANSMISSION_NETWORKS_VERSION_H
#define TRANSMISSION_NETWORKS_VERSION_H

#include <string>

namespace transmission_nets::core::version {

// Version information
// These are set by CMake during build configuration
#ifndef TRANSMISSION_NETWORKS_VERSION_MAJOR
#define TRANSMISSION_NETWORKS_VERSION_MAJOR 1
#endif

#ifndef TRANSMISSION_NETWORKS_VERSION_MINOR
#define TRANSMISSION_NETWORKS_VERSION_MINOR 0
#endif

#ifndef TRANSMISSION_NETWORKS_VERSION_PATCH
#define TRANSMISSION_NETWORKS_VERSION_PATCH 0
#endif

// Version string
inline std::string getVersionString() {
    return std::to_string(TRANSMISSION_NETWORKS_VERSION_MAJOR) + "." +
           std::to_string(TRANSMISSION_NETWORKS_VERSION_MINOR) + "." +
           std::to_string(TRANSMISSION_NETWORKS_VERSION_PATCH);
}

// Version components
inline int getMajorVersion() { return TRANSMISSION_NETWORKS_VERSION_MAJOR; }
inline int getMinorVersion() { return TRANSMISSION_NETWORKS_VERSION_MINOR; }
inline int getPatchVersion() { return TRANSMISSION_NETWORKS_VERSION_PATCH; }

} // namespace transmission_nets::core::version

#endif // TRANSMISSION_NETWORKS_VERSION_H
