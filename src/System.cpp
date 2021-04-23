#include "System.h"
#include "utils.h"
#include <cmath>
#include <tuple>


System::System(Jet *newjet,
               Vector3d &newpoint_in,
               Vector3d &newray_direction,
               double newnu) {
    jet = newjet;
    point_in = newpoint_in;
    ray_direction = newray_direction;
    nu = newnu;
}


void Tau::operator()(const double &x, double &dxdt, const double t) {
    Vector3d point = point_in - t * ray_direction;
    dxdt = jet->getKI(point, ray_direction, nu);
}


void TauFR::operator()(const double &x, double &dxdt, double t) {
    Vector3d point = point_in + t * ray_direction;
    dxdt = jet->getKF(point, ray_direction, nu);
}


void I::operator()(const double &x, double &dxdt, const double t) {
    Vector3d point = point_in + t * ray_direction;

    // If negative Stokes I resulted from previous step, remember its value to compensate for
    double compensate_negative_I = 0.0;
    if(x < 0) {
        compensate_negative_I = -x;
    }

    double k_i, eta_i;
    std::tie(k_i, eta_i) = jet->get_stokes_I_transport_coefficients(point, ray_direction, nu);

    dxdt = eta_i - k_i*(x+compensate_negative_I);

    // FIXME: This is possibly wrong
    // This adds to previous step Stokes I, so add value that compensating negative Stokes I from previous step
    dxdt += compensate_negative_I;
}


void Speed::operator()(const double &x, double &dxdt, const double t) {
    Vector3d point = point_in + t * ray_direction;

    // If negative Stokes I resulted from previous step, remember its value to compensate for
    double compensate_negative_I = 0.0;
    if(x < 0) {
        compensate_negative_I = -x;
    }

    double k_i, eta_i;
    std::tie(k_i, eta_i) = jet->get_stokes_I_transport_coefficients(point, ray_direction, nu);

    dxdt = eta_i - k_i*(x+compensate_negative_I);

    // This adds to previous step Stokes I, so add value that compensating negative Stokes I from previous step
    dxdt += compensate_negative_I;

    // Calculate apparent speed
    double psi = jet->getPsi(point);
    Vector3d beta = jet->getV(point, psi)/c;
    // LOS angle is angle between the bulk motion speed and the local ``ray_direction`` (that goes from far side of the
    // source to the observer).
    double cosLOS = beta.dot(ray_direction)/beta.norm();
    double sinLOS = sqrt(1 - cosLOS*cosLOS);
//    std::cout << "cosLOS = " << cosLOS << ", sinLOS = " << sinLOS << "\n";
    double beta_app = beta.norm()*sinLOS / (1. - beta.norm()*cosLOS);
//    std::cout << "beta_app = " << beta_app << "\n";
    if(isnan(beta_app)) {
        beta_app = 0.0;
//        std::cout << "beta_app = NAN!" << "\n";
//        std::cout << "sinLOS = " << sinLOS << "\n";
//        std::cout << "beta = " << beta << "\n";
//        std::cout << "Psi = " << psi << "\n";
    }

    dxdt *= beta_app;
}


FullStokes::FullStokes(Jet *newjet, Vector3d &newpoint_in,
                       Vector3d &newray_direction, double newnu) {
	jet = newjet;
	point_in = newpoint_in;
	ray_direction = newray_direction;
	nu = newnu;
}

void FullStokes::operator()(const state_type &x, state_type &dxdt,
                            const double t) {
    Vector3d point = point_in + t * ray_direction;

	// If negative Stokes I resulted from previous step, remember its value to compensate for
	double compensate_negative_I = 0.0;
	if(x[0] < 0) {
	     compensate_negative_I = -x[0];
	}



	// TODO: Do I need this block of checking non-neglible B-field?
	double psi = jet->getPsi(point);
	// Find value of B in lab frame to decide about making transport or not
	// FIXME: This is not always plasma-frame B!
	auto b_in_plasma_frame = jet->getB(point, psi);
	auto bvalues_in_plasma_frame = b_in_plasma_frame.norm();



	// Non-zero B => make some transport!
	if(bvalues_in_plasma_frame > eps_B) {

        Vector3d B_hat = jet->getBhat(point, psi);
        Vector3d beta = jet->getV(point, psi)/c;
        // Unit vector of electric vector of wave (EVPA) in lab frame.
        // FIXME: If beta = 0 => e must be perpendicular to B_hat --- fixed inside ``e_hat`` (in function ``utils.q()``)
        Vector3d e = e_hat(B_hat, beta, ray_direction);

        // Reference direction l (on sky along the jet axis). When chi = 0.5*atan2(U,Q) = 0 than EVPA should be along the
        // jet.
        Vector3d l_y(0.0, 1.0, 0.0);
        auto l = ray_direction.cross(l_y);

        // Need inverse rotation (on -2chi) from ``observers system`` to ``emitting system`` (see Porth's PhD)
        auto cos_chi = l.dot(e);
        auto sin_chi = ray_direction.dot(l.cross(e));
        auto sin_2chi = 2.0*sin_chi*cos_chi;
        auto cos_2chi = cos_chi*cos_chi - sin_chi*sin_chi;
        // ``x_eb`` is State in ``emitting system``, ``x`` - in `observers system``
        auto x_eb = x;
        x_eb[1] = x[1]*cos_2chi - x[2]*sin_2chi;
        x_eb[2] = x[1]*sin_2chi + x[2]*cos_2chi;
        // Previous step compensated for negative Stokes I flux - to calculate change of state dxdt
        x_eb[0] += compensate_negative_I;

        // Change of state in ``emitting system``
        auto dxdt_eb = dxdt;

        double k_i;
        double k_q;
        double k_u;
        double k_v;
        double eta_i;
        double eta_q;
        double eta_u;
        double eta_v;
        double k_F;
        double k_C;
        double h_Q;

        std::tie(k_i, k_q, k_u, k_v, eta_i, eta_q, eta_u, eta_v, k_F, k_C, h_Q) = jet->get_transport_coefficients(point, ray_direction, nu);

        dxdt_eb[0] = eta_i - k_i*x_eb[0] - k_q*x_eb[1] - k_u*x_eb[2] - k_v*x_eb[3];
        dxdt_eb[1] = eta_q - k_i*x_eb[1] - k_q*x_eb[0] - k_F*x_eb[2] - h_Q*x_eb[3];
        dxdt_eb[2] = eta_u - k_i*x_eb[2] - k_u*x_eb[0] + k_F*x_eb[1] - k_C*x_eb[3];
        dxdt_eb[3] = eta_v - k_i*x_eb[3] - k_v*x_eb[0] + h_Q*x_eb[1] + k_C*x_eb[2];

        // Need rotation (on 2chi) from ``emitting system` to ``observers system``
        // Change of state in ``observers system``
        dxdt = dxdt_eb;

        // FIXME: This is possibly wrong
        // This adds to previous step Stokes I, so add value that compensating negative Stokes I from previous step
        // FIXME: I must add to x[0] (however it is const)!!!
        dxdt[0] += compensate_negative_I;
        dxdt[1] = dxdt_eb[1]*cos_2chi + dxdt_eb[2]*sin_2chi;
        dxdt[2] = -dxdt_eb[1]*sin_2chi + dxdt_eb[2]*cos_2chi;

        if (isnan(dxdt.at(0)) || isnan(dxdt.at(1)) || isnan(dxdt.at(2)) || isnan(dxdt.at(3))) {
            std::cout << "Problems in transfer!" << std::endl;
        }
    }

	// Zero B => no changes in radiation at all!
	else {
        dxdt[0] = 0.0;
        dxdt[1] = 0.0;
        dxdt[2] = 0.0;
        dxdt[3] = 0.0;
	}

}

bool check_opt_depth(double tau_max, const double &x) {
    return x >= tau_max;
}





//class push_back_state_and_time
//{
//    public:
//        push_back_state_and_time(std::vector<state_type> &states,
//                                 std::vector<double> &times,
//                                 std::vector<double> &P_0,
//                                 Jet* new_jet) :
//                counter(0),
//                m_states(states),
//                m_times(times),
//                m_P_cur(P_0)
//        {
//            jet = new_jet;
//        }
//        void operator()(const state_type &x, double t)
//        {
//            if(std::isnan(x[0])) {
//                throw JetBorderException();
//            }
//            //m_states.clear();
//            //m_times.clear();
//            m_states.push_back(x);
//            m_times.push_back(t);
//
//            double p_cur = solver->get_P0(x, t);
//            m_P_cur.push_back(p_cur);
//
//        }
//    private:
//        std::vector<state_type>& m_states;
//        std::vector<double>& m_times;
//        std::vector<double>& m_P_cur;
//        Jet* jet;
//        size_t counter;
//};
