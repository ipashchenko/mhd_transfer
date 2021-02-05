#include <cmath>
#include <Eigen/Eigen>
#include <utils.h>
#include "VField.h"

using Eigen::Vector3d;

SimulationVField::SimulationVField(Delaunay_triangulation *tr) : interp_(tr, 1.0) {}

Vector3d SimulationVField::vf(const Vector3d &point) const {
    double gamma = interp_.interpolated_value(point);
    return Vector3d{0, 0, c * sqrt(1. - 1. / (gamma * gamma))};
};
