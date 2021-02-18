#include <cmath>
#include <Eigen/Eigen>
#include <BField.h>
#include <utils.h>
#include <boost/math/special_functions/bessel.hpp>

using Eigen::Vector3d;


VectorBField::VectorBField(bool in_plasma_frame, double tangled_fraction, Geometry* geometry_out, Geometry* geometry_in) :
        in_plasma_frame_(in_plasma_frame),
        tangled_fraction_(tangled_fraction) {
    geometry_in_ = geometry_in;
    geometry_out_ = geometry_out;
}

Vector3d VectorBField::bf(const Vector3d &point, double psi) const {
    double x = point[0];
    double y = point[1];
    double r_point = sqrt(x*x + y*y);

    if(geometry_out_) {
        // Find radius of outer surface at given point
        double r_border_out = geometry_out_->radius_at_given_distance(point);
        if (r_point > r_border_out) {
            return {0.0, 0.0, 0.0};
        }
    }
    if(geometry_in_) {
        // Find radius of inner surface at given point
        double r_border_in = geometry_in_->radius_at_given_distance(point);
        if (r_point < r_border_in) {
            return {0.0, 0.0, 0.0};
        }
    }
    return _bf(point, psi);
}

Vector3d VectorBField::bf_plasma_frame(const Vector3d &point, double psi, Vector3d &v) const {
    Vector3d b = bf(point, psi);
    if (in_plasma_frame_) {
        return b;
    } else {
        return get_B_prime(b, v);
    }
}

double VectorBField::bf_tangled_plasma_frame(const Vector3d &point, double psi, Vector3d &v) const {
    Vector3d b = bf(point, psi);
    if (in_plasma_frame_) {
        return tangled_fraction_*b.norm();
    } else {
        return tangled_fraction_*get_B_prime(b, v).norm();
    }
}

Vector3d VectorBField::bhat_lab_frame(const Vector3d &point, double psi, Vector3d &v) const {
    Vector3d b = bf(point, psi);
    if (in_plasma_frame_) {
        auto b_hat_prime = b.normalized();
        Vector3d beta = v/c;
        return get_B_hat(b_hat_prime, beta);
    } else {
        return b.normalized();
    }
}

SimulationBField::SimulationBField(Delaunay_triangulation *tr_p, Delaunay_triangulation *tr_fi, bool in_plasma_frame, double tangled_fraction) :
        VectorBField(in_plasma_frame, tangled_fraction, nullptr, nullptr), interp_p_(tr_p), interp_fi_(tr_fi) {}

Vector3d SimulationBField::_bf(const Vector3d &point, double psi) const {
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double phi = atan2(y, x);
    // As atan2 returns [-pi, pi], put this to [0, 2pi]
    if(phi < 0){
        phi += 2.0*M_PI;
    }
    double interpolated_value_p = interp_p_.interpolated_value({psi, z/pc});
    double interpolated_value_fi = interp_fi_.interpolated_value({psi, z/pc});
    return Vector3d{-sin(phi) * interpolated_value_fi, cos(phi) * interpolated_value_fi, interpolated_value_p};
}