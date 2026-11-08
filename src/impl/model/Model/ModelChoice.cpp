//
// The Run recipe's model choices, parsed and in canonical form. See ModelChoice.h.
//

#include "ModelChoice.h"

#include <fmt/core.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace transmission_nets::impl::Model {
    namespace {
        // The number `text` spells, all of it.
        double wholeNumber(const std::string& text) {
            std::size_t used = 0;
            const double x   = std::stod(text, &used);
            if (used != text.size()) {
                throw std::invalid_argument(text);
            }
            return x;
        }
    }// namespace

    double LeakPrior::logDensity(const double leak) const {
        if (kind == Kind::LogUniform) {
            return leak >= lo and leak <= hi ? -std::log(leak) - std::log(std::log(hi / lo)) : -std::numeric_limits<double>::infinity();
        }
        return leak > 0 and leak < EXPONENTIAL_UPPER ? std::log(rate) - rate * leak : -std::numeric_limits<double>::infinity();
    }

    double LeakPrior::startAt(const double leak) const {
        if (holds(leak)) {
            return leak;
        }
        return std::clamp(leak, lower() * 1.0001 + 1e-12, upper() * 0.9999);
    }

    LeakPrior parseLeakPrior(const std::string& description) {
        const auto colon       = description.find(':');
        const std::string kind = description.substr(0, colon);
        const std::string args = colon == std::string::npos ? "" : description.substr(colon + 1);
        try {
            if (kind == "loguniform") {
                LeakPrior prior{.kind = LeakPrior::Kind::LogUniform};
                if (not args.empty()) {
                    const auto comma = args.find(',');
                    if (comma == std::string::npos) {
                        throw std::invalid_argument(args);
                    }
                    prior.lo = wholeNumber(args.substr(0, comma));
                    prior.hi = wholeNumber(args.substr(comma + 1));
                }
                if (prior.lo > 0 and prior.lo < prior.hi and prior.hi < 1) {
                    return prior;
                }
            } else if (kind == "exponential") {
                LeakPrior prior{.kind = LeakPrior::Kind::Exponential};
                if (not args.empty()) {
                    prior.rate = wholeNumber(args);
                }
                if (prior.rate > 0) {
                    return prior;
                }
            }
        } catch (const std::exception&) {
        }
        throw std::invalid_argument(fmt::format("a leak prior is loguniform:LO,HI with 0 < LO < HI < 1, or exponential:RATE with RATE > 0 (got \"{}\")", description));
    }

    std::string canonicalLeakPrior(const std::string& description) {
        try {
            const auto prior = parseLeakPrior(description);
            return prior.kind == LeakPrior::Kind::LogUniform ? fmt::format("loguniform:{},{}", prior.lo, prior.hi) : fmt::format("exponential:{}", prior.rate);
        } catch (const std::invalid_argument&) {
            return description;
        }
    }

    ErrorModel parseErrorModel(const std::string& name) {
        if (name.empty() or name == "infection") {
            return ErrorModel::Infection;
        }
        if (name == "shared") {
            return ErrorModel::Shared;
        }
        if (name == "hierarchical") {
            return ErrorModel::Hierarchical;
        }
        throw std::invalid_argument(fmt::format("an error model is infection, shared or hierarchical (got \"{}\")", name));
    }

    std::string canonicalErrorModel(const std::string& name) {
        try {
            return parseErrorModel(name) == ErrorModel::Infection ? std::string{} : name;
        } catch (const std::invalid_argument&) {
            return name;
        }
    }
}// namespace transmission_nets::impl::Model
