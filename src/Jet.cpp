#include <tuple>
#include "Jet.h"
#include "MyExceptions.h"

using Eigen::Vector3d;


Jet::Jet(BaseGeometry *newgeo, SimulationInterpolater *newPsiInterpolater, VField *newvfield, VectorBField* newbField, NField* newnField) {
    geometry_ = newgeo;
    PsiInterpolater_ = newPsiInterpolater;
    vfield_ = newvfield;
    bfield_ = newbField;
    nfield_ = newnField;
}


std::tuple<double, double, double, double, double, double, double, double, double, double, double> Jet::get_transport_coefficients(Vector3d &point, Vector3d &n_los, double nu) {
    // Example for k_I (Lyutikov et al. 2005):
    // First, comoving frame ``k_i_prime`` (in the rest frame of the emission element) is connected to this ``k_i`` as
    // ``k_i = k_i_prime / D``. Second, in ``k_i_prime`` we need all quantities in comoving frame (primed) in terms of
    // lab frame:
    // b_prime = f(b, v)
    // n_los_prime = f(n_los, v)
    // nu_prime = f(nu, n_los, v) = nu/getD
    // n_prime = f(n, v) = n/Gamma

    // Get Psi given (r, z)
    double psi = getPsi(point);

    Vector3d v = getV(point, psi);
    auto gamma = getG(v);

    auto b_prime = bfield_->bf_plasma_frame(point, psi, v);

    //if(b_prime.norm() < eps_B || isnan(b_prime.norm())) {
    //    return std::make_tuple(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    //}
    // NaN means something is wrong!
    if(b_prime.norm() < eps_B) {
        return std::make_tuple(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }

    auto b_prime_tangled = b_prime.norm()*bfield_->get_tangled_fraction();

    auto D = getD(n_los, v);
    auto nu_prime = nu/D;
    auto n_los_prime = get_n_los_prime(n_los, v);

    // Now calculate all coefficients
    double k_i_prime = 0.0;
    double k_q_prime = 0.0;
    double k_u_prime = 0.0;
    double k_v_prime = 0.0;
    double eta_i_prime = 0.0;
    double eta_q_prime = 0.0;
    double eta_u_prime = 0.0;
    double eta_v_prime = 0.0;
    double k_F_prime = 0.0;
    double k_C_prime = 0.0;
    double h_Q_prime = 0.0;
    double n_prime = nfield_->nf_plasma_frame(point, psi, gamma);
    k_i_prime += nfield_->k_i(b_prime, n_los_prime, nu_prime, n_prime);
    k_i_prime += nfield_->k_i(b_prime_tangled, n_los_prime, nu_prime, n_prime);
    k_q_prime += nfield_->k_q(b_prime, n_los_prime, nu_prime, n_prime);
    k_u_prime += nfield_->k_u(b_prime, n_los_prime, nu_prime, n_prime);
    k_v_prime += nfield_->k_v(b_prime, n_los_prime, nu_prime, n_prime);
    eta_i_prime += nfield_->eta_i(b_prime, n_los_prime, nu_prime, n_prime);
    eta_i_prime += nfield_->eta_i(b_prime_tangled, n_los_prime, nu_prime, n_prime);
    eta_q_prime += nfield_->eta_q(b_prime, n_los_prime, nu_prime, n_prime);
    eta_u_prime += nfield_->eta_u(b_prime, n_los_prime, nu_prime, n_prime);
    eta_v_prime += nfield_->eta_v(b_prime, n_los_prime, nu_prime, n_prime);
    k_F_prime += nfield_->k_F(b_prime, n_los_prime, nu_prime, n_prime);
    k_C_prime += nfield_->k_C(b_prime, n_los_prime, nu_prime, n_prime);
    h_Q_prime += nfield_->h_Q(b_prime, n_los_prime, nu_prime, n_prime);

    //std::cout << "k_F = " << k_F_prime/D << std::endl;
    auto result = std::make_tuple(k_i_prime/D, k_q_prime/D, k_u_prime/D, k_v_prime/D,
                                  eta_i_prime*D*D, eta_q_prime*D*D, eta_u_prime*D*D, eta_v_prime*D*D,
                                  k_F_prime/D, k_C_prime/D, h_Q_prime/D);
    if(isnan(k_i_prime/D)) {
        std::cout << "NaN in k_I!" << std::endl;
    }
    return result;
}


std::tuple<double, double> Jet::get_stokes_I_transport_coefficients(Vector3d &point, Vector3d &n_los, double nu) {
    // Example for k_I (Lyutikov et al. 2005):
    // First, comoving frame ``k_i_prime`` (in the rest frame of the emission element) is connected to this ``k_i`` as
    // ``k_i = k_i_prime / D``. Second, in ``k_i_prime`` we need all quantities in comoving frame (primed) in terms of
    // lab frame:
    // b_prime = f(b, v)
    // n_los_prime = f(n_los, v)
    // nu_prime = f(nu, n_los, v) = nu/getD
    // n_prime = f(n, v) = n/Gamma

    // Get Psi given (r, z)
    double psi = getPsi(point);

    Vector3d v = getV(point, psi);
    auto gamma = getG(v);

    auto b_prime = bfield_->bf_plasma_frame(point, psi, v);

    if(b_prime.norm() < eps_B) {
        return std::make_tuple(0.0, 0.0);
    }

    auto b_prime_tangled = b_prime.norm()*bfield_->get_tangled_fraction();

    auto D = getD(n_los, v);
    auto nu_prime = nu/D;
    auto n_los_prime = get_n_los_prime(n_los, v);

    // Now calculate all coefficients
    double k_i_prime = 0.0;
    double eta_i_prime = 0.0;
    double n_prime = nfield_->nf_plasma_frame(point, psi, gamma);
    k_i_prime += nfield_->k_i(b_prime, n_los_prime, nu_prime, n_prime);
    k_i_prime += nfield_->k_i(b_prime_tangled, n_los_prime, nu_prime, n_prime);
    eta_i_prime += nfield_->eta_i(b_prime, n_los_prime, nu_prime, n_prime);
    eta_i_prime += nfield_->eta_i(b_prime_tangled, n_los_prime, nu_prime, n_prime);

    std::cout << "=== k, eta = " << k_i_prime/D << ", " << eta_i_prime*D*D << "\n";

    return std::make_tuple(k_i_prime/D, eta_i_prime*D*D);
//    return std::make_tuple(1E+23*k_i_prime/D, 1e+23*eta_i_prime*D*D);
}




// This is k_i in lab frame that could be integrated along LOS.
double Jet::getKI(Vector3d &point, Vector3d &n_los, double nu) {
    // First, comoving frame ``k_i_prime`` (in the rest frame of the emission
    // element) is connected to this ``k_i`` as ``k_i = k_i_prime / D``.
    // Second, in ``k_i_prime`` we need all quantities in comoving frame
    // (primed) in terms of lab frame:
    // b_prime = f(b, v)
    // n_los_prime = f(n_los, v)
    // nu_prime = f(nu, n_los, v) = nu/getD
    // n_prime = f(n, v) = n/Gamma

    // Get Psi given (r, z)
    double psi = getPsi(point);

    Vector3d v = getV(point, psi);
    auto gamma = getG(v);

    auto b_prime = bfield_->bf_plasma_frame(point, psi, v);

    //if(b_prime.norm() < eps_B || isnan(b_prime.norm())) {
    //    return 0.0;
    //}
    // NaN means something is wrong!
    if(b_prime.norm() < eps_B) {
        return 0.0;
    }

    auto b_prime_tangled = b_prime.norm()*bfield_->get_tangled_fraction();

    auto D = getD(n_los, v);
    auto nu_prime = nu/D;
    auto n_los_prime = get_n_los_prime(n_los, v);

    double k_i_prime = 0.0;
    double n_prime = nfield_->nf_plasma_frame(point, psi, gamma);
    k_i_prime += nfield_->k_i(b_prime, n_los_prime, nu_prime, n_prime);
    k_i_prime += nfield_->k_i(b_prime_tangled, n_los_prime, nu_prime, n_prime);

    auto result = k_i_prime/D;
    if(isnan(result)) {
        std::cout << "NaN in k_I!" << std::endl;
    }
    std::cout << "k_I = " << result << "\n";
    return result;
}

// This is eta_i in lab frame that could be integrated along LOS.
double Jet::getEtaI(Vector3d &point, Vector3d &n_los, double nu) {
    // First, comoving frame ``eta_i_prime`` (in the rest frame of the emission
    // element) is connected to this ``eta_i`` as ``eta_i = D^2 * eta_i_prime``.
    // Second, in ``eta_i_prime`` we need all quantities in comoving frame
    // (primed) in terms of lab frame:
    // b_prime = f(b, v)
    // n_los_prime = f(n_los, v)
    // nu_prime = f(nu, n_los, v) = nu/getD
    // n_prime = f(n, v) = n/Gamma

    // Get Psi given (r, z)
    double psi = getPsi(point);

    Vector3d v = getV(point, psi);
    auto gamma = getG(v);

    auto b_prime = bfield_->bf_plasma_frame(point, psi, v);

    if(b_prime.norm() < eps_B || isnan(b_prime.norm())) {
        return 0.0;
    }

    auto b_prime_tangled = b_prime.norm()*bfield_->get_tangled_fraction();

    auto D = getD(n_los, v);
    auto nu_prime = nu/D;
    auto n_los_prime = get_n_los_prime(n_los, v);

    double eta_i_prime = 0.0;
    double n_prime = nfield_->nf_plasma_frame(point, psi, gamma);
    eta_i_prime += nfield_->eta_i(b_prime, n_los_prime, nu_prime, n_prime);
    eta_i_prime += nfield_->eta_i(b_prime_tangled, n_los_prime, nu_prime, n_prime);

    auto result = eta_i_prime*D*D;
    if(isnan(result)) {
        std::cout << "NaN in eta_I!" << std::endl;
    }
    return result;
}


double Jet::getKF(Vector3d &point, Vector3d &n_los, double nu) {

    // Get Psi given (r, z)
    double psi = getPsi(point);

    Vector3d v = getV(point, psi);
    auto gamma = getG(v);

    auto b_prime = bfield_->bf_plasma_frame(point, psi, v);

    if(b_prime.norm() < eps_B || isnan(b_prime.norm())) {
        return 0.0;
    }

    auto D = getD(n_los, v);
    auto nu_prime = nu/D;
    auto n_los_prime = get_n_los_prime(n_los, v);

    double k_F_prime = 0.0;
    double n_prime = nfield_->nf_plasma_frame(point, psi, gamma);
    k_F_prime += nfield_->k_F(b_prime, n_los_prime, nu_prime, n_prime);

    auto result = k_F_prime/D;
    //std::cout << "k_F = " << result << std::endl;
    if(isnan(result)) {
        std::cout << "NaN in k_F!" << std::endl;
    }
    return result;
}


std::list<Intersection> Jet::hit(Ray &ray) {
    return geometry_->hit(ray);
}


double Jet::getPsi(const Vector3d &point) {
    double x = point[0]/pc;
    double y = point[1]/pc;
    double z = point[2]/pc;
    return PsiInterpolater_->interpolated_value({hypot(x, y), z});
}

Vector3d Jet::getV(const Vector3d &point, double psi) {
    auto v = vfield_->vf(point, psi);
    if(v.norm() > c) {
        std::cout << "Speed > c!!!";
        throw PhysicalException("Speed");
    }
    return v;
}

const Vector3d Jet::getB(const Vector3d &point, double psi) {
    return bfield_->bf(point, psi);
}

const Vector3d Jet::getBhat(const Vector3d& point, double psi) {
    auto v = getV(point, psi);
    return bfield_->bhat_lab_frame(point, psi, v);
}
