//
// Created by Maxwell Murphy on 1/23/20.
//

#include "LogLikelihood.h"



LogLikelihood::LogLikelihood(float value, std::string id) : value_(value), id_(std::move(id)) {}

LogLikelihood::LogLikelihood(std::string id) : value_(0), id_(std::move(id)) {}
