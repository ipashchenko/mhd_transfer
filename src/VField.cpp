#include <cmath>
#include <Eigen/Eigen>
#include <utils.h>
#include "VField.h"

using Eigen::Vector3d;

SimulationVField::SimulationVField(Delaunay_triangulation *tr_Gamma, Delaunay_triangulation *tr_beta_phi) :
    interp_Gamma_(tr_Gamma, 1.0), interp_beta_phi_(tr_beta_phi, 0.0) {}

Vector3d SimulationVField::vf(const Vector3d &point) const {

    double x = point[0];
    double y = point[1];
    double z = point[2];
    double phi = atan2(y, x);

    double gamma = interp_Gamma_.interpolated_value(point);
    double v_phi = c*interp_beta_phi_.interpolated_value(point);

    if(z > 0) {
        return {-sin(phi)*v_phi/gamma, cos(phi)*v_phi/gamma, c*sqrt(1. - 1./(gamma*gamma))};
        //return {sin(phi)*v_phi, -cos(phi)*v_phi, c*sqrt(1. - 1./(gamma*gamma))};
    } else {
        return {sin(phi)*v_phi/gamma, -cos(phi)*v_phi/gamma, -c*sqrt(1. - 1./(gamma*gamma))};
    }
    return Vector3d{0, 0, c * sqrt(1. - 1. / (gamma * gamma))};
};
