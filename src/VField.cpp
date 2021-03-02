#include <cmath>
#include <Eigen/Eigen>
#include <utils.h>
#include "VField.h"

using Eigen::Vector3d;

SimulationVField::SimulationVField(Delaunay_triangulation *tr_psi, Delaunay_triangulation *tr_Gamma, Delaunay_triangulation *tr_beta_phi) :
        interp_psi_(tr_psi, 1.0), interp_Gamma_(tr_Gamma, 1.0), interp_beta_phi_(tr_beta_phi, 0.0) {}

Vector3d SimulationVField::vf(const Vector3d &point, double psi) const {

    double x = point[0];
    double y = point[1];
    double z = point[2];
    double phi = atan2(y, x);

    // As atan2 returns [-pi, pi], put this to [0, 2pi]
    if(phi < 0){
        phi += 2.0*M_PI;
    }

    // It is total speed!
    double gamma = interp_Gamma_.interpolated_value({psi, z/pc});
    double v_tot = c*sqrt(gamma*gamma - 1.0)/gamma;
    double v_phi = c*interp_beta_phi_.interpolated_value({psi, z/pc});
    double v_pol = sqrt(v_tot*v_tot - v_phi*v_phi);

    // Find direction of the poloidal component
    Vector2d gradPsi = interp_psi_.gradient({hypot(x, y)/pc, z/pc});
    gradPsi.normalize();
    // sin & cos of angle between the poloidal component and z-axis
    double sinz = gradPsi[0];
    double cosz = gradPsi[1];
    Vector3d V_p;
    if(z > 0 ){
        V_p = {v_pol*sinz*cos(phi), v_pol*sinz*sin(phi), v_pol*cosz};
    } else {
        V_p = {v_pol*sinz*cos(phi), v_pol*sinz*sin(phi), -v_pol*cosz};
    }

    Vector3d V_phi = {-sin(phi)*v_phi, cos(phi)*v_phi, 0.0};

    return V_p + V_phi;
//    if(z > 0) {
//        return {-sin(phi)*v_phi, cos(phi)*v_phi, v_pol};
//    } else {
//        return {sin(phi)*v_phi, -cos(phi)*v_phi, -v_pol};
//    }
};
