#pragma once
#include "Types.hpp"
#include "basis.h"
#include <Eigen/Core>

namespace SIBSplines
{

std::vector<int> feasible_control_point_of_given_parameter(const double para,
                                                           const std::vector<double> &U,
                                                           const int degree, const double per);
}