#include <cmath>
#include <Eigen/Eigen>
#include <utils.h>
#include "VField.h"

using Eigen::Vector3d;

SimulationVField::SimulationVField(Delaunay_triangulation *tr_Gamma, Delaunay_triangulation *tr_beta_phi) :
    interp_Gamma_(tr_Gamma, 1.0), interp_beta_phi_(tr_beta_phi, 0.0) {}

Vector3d SimulationVField::vf(const Vector3d &point, double psi) const {

    double x = point[0];
    double y = point[1];
    double z = point[2];
    double phi = atan2(y, x);
    // As atan2 returns [-pi, pi], put this to [0, 2pi]
    if(phi < 0){
        phi += 2.0*M_PI;
    }

    double gamma = interp_Gamma_.interpolated_value({psi, z/pc});
    double v_phi = c*interp_beta_phi_.interpolated_value({psi, z/pc});

    if(z > 0) {
        return {-sin(phi)*v_phi, cos(phi)*v_phi, c*sqrt(1. - 1./(gamma*gamma))};
    } else {
        return {sin(phi)*v_phi, -cos(phi)*v_phi, -c*sqrt(1. - 1./(gamma*gamma))};
    }
};
