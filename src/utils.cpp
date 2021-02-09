#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include "utils.h"
#include "MyExceptions.h"


double nu_p(double n) {return sqrt(n*q_e*q_e / (pi*m_e));}

double nu_b(Vector3d &b, Vector3d &n_los) {
    return q_e*(n_los.cross(b)).norm()/(2.*pi*m_e*c);
}

double nu_b(double b) {
	return q_e*b/(2.*pi*m_e*c);
}

double nu_min(double gamma_min, Vector3d &b, Vector3d &n_los) {
    return gamma_min*gamma_min*nu_b(b, n_los);
}

double nu_min(double gamma_min, double b) {
    return gamma_min*gamma_min*nu_b(b);
}

double sin_theta(Vector3d &b, Vector3d &n_los) {
		return n_los.cross(b).norm()/b.norm();
};

double nu_b_value(Vector3d &b) {
	return q_e*b.norm()/(2.*pi*m_e*c);
}

double k_0(Vector3d &b, Vector3d &n_los, double nu, double n) {
    return pi*nu_p(n)*nu_p(n)*nu_b(b, n_los)/(c*nu*nu);
}

double k_0(double b, Vector3d &n_los, double nu, double n) {
	return pi*nu_p(n)*nu_p(n)*nu_b(b)/(c*nu*nu);
}

double k_0_value(Vector3d &b, double nu, double n) {
	return pi*nu_p(n)*nu_p(n)*nu_b_value(b)/(c*nu*nu);
}

double k_i(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
    double factor = (pow(3., (s+1.)/2.)/4.)*tgamma(s/4.+11./6.)*tgamma(s/4.+1./6.);
    return k_0(b, n_los, nu, n) * pow(nu_b(b, n_los)/nu, s/2.) * factor;
}

double k_i(double b, Vector3d &n_los, double nu, double n, double s) {
	double factor = (pow(3., (s+1.)/2.)/4.)*tgamma(s/4.+11./6.)*tgamma(s/4.+1./6.);
	double rnd_factor = sqrt(pi/4.)*tgamma((6.+s)/4.)/tgamma((8.+s)/4.);
	factor = factor*rnd_factor;
	return k_0(b, n_los, nu, n) * pow(nu_b(b)/nu, s/2.) * factor;
}

// Alternative formulation
double k_i_(double b, Vector3d &n_los, double nu, double n, double s) {
	double alpha = (s-1.)/2.;
	double factor = sqrt(pi/4.)*tgamma((7.+2.*alpha)/4.)/tgamma((9.+2.*alpha)/4.);
	factor = factor*pow(m_e*c*c, 2.*alpha)*(sqrt(3.)*pow(q_e, 3.)/(8.*pi*m_e));
	factor = factor*pow((3.*q_e)/(2.*pi*pow(m_e, 3.)*pow(c, 5.)), alpha+0.5);
	factor = factor*tgamma((6.*alpha+5.)/12.)*tgamma((6.*alpha+25.)/12.);
	return factor*pow(b, alpha+1.5)*n*pow(nu, -alpha-2.5);
}

double k_q(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
	return (s+2.)/(s+10./3)*k_i(b, n_los, nu, n, s);
}

double k_u(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
	return 0;
}

double k_v(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
	double cos_theta = b.dot(n_los)/b.norm();
	double factor = pow(3., s/2.)*(s+3.)*(s+2.)/(4.*(s+1.))*tgamma(s/4.+11./12.)*tgamma(s/4.+7./12.);
	return -k_0_value(b, nu, n)*cos_theta*pow(nu_b(b, n_los)/nu, (s+1.)/2.) * factor;
}

double k_F_c(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
	double cos_theta = b.dot(n_los)/b.norm();
	return 2.*pi*nu_p(n)*nu_p(n)*nu_b_value(b)*cos_theta/(c*nu*nu);
}

double k_C_c(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
	return -pi*nu_p(n)*nu_p(n)*nu_b(b, n_los)*nu_b(b, n_los)/(c*nu*nu*nu);
}

double k_F_r(Vector3d &b, Vector3d &n_los, double nu, double n, double gamma_min, double s) {
	return (s+2.)*log(gamma_min)/((s+1.)*pow(gamma_min, s+1.))*k_F_c(b, n_los, nu, n, s);
}

double k_C_r(Vector3d &b, Vector3d &n_los, double nu, double n, double gamma_min, double s) {
	return (2./(s-2.))*(pow(gamma_min, 2.-s) - pow(nu_b(b, n_los)/nu, (s-2.)/2.))*k_C_c(b, n_los, nu, n, s);
}

double h_Q(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
	return 0;
}

double eta_0(Vector3d &b, Vector3d &n_los, double n) {
    return pi*nu_p(n)*nu_p(n)*nu_b(b, n_los)*m_e/c;
}

double eta_0(double b, Vector3d &n_los, double n) {
	return pi*nu_p(n)*nu_p(n)*nu_b(b)*m_e/c;
}

double eta_0_value(Vector3d &b, double n) {
	return pi*nu_p(n)*nu_p(n)*nu_b_value(b)*m_e/c;
}

double eta_i(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
    double factor = pow(3., s/2.)/(2.*(s+1))*tgamma(s/4.+19./12.)*tgamma(s/4.-1./12.);
    return eta_0(b, n_los, n) * pow(nu_b(b, n_los)/nu, (s-1.)/2.) * factor;
}

double eta_i(double b, Vector3d &n_los, double nu, double n, double s) {
	double factor = pow(3., s/2.)/(2.*(s+1))*tgamma(s/4.+19./12.)*tgamma(s/4.-1./12.);
	double rnd_factor = sqrt(pi/4.)*tgamma((5.+s)/4.)/tgamma((7.+s)/4.);
	factor = factor*rnd_factor;
	return eta_0(b, n_los, n) * pow(nu_b(b)/nu, (s-1.)/2.) * factor;
}

double eta_q(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
	return (s+1.0)/(s+7./3.)*eta_i(b, n_los, nu, n, s);
}

double eta_u(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
	return 0;
}

double eta_v(Vector3d &b, Vector3d &n_los, double nu, double n, double s) {
	double cos_theta = b.dot(n_los)/b.norm();
	double factor = pow(3., (s-1.)/2.)*(s+2.)/(2.*s)*tgamma(s/4.+2./3.)*tgamma(s/4.+1./3.);
	return -eta_0_value(b, n)*cos_theta*pow(nu_b(b, n_los)/nu, s/2.) * factor;
}

double getG(Vector3d &v) {
    Vector3d beta = v/c;
    return sqrt(1./(1.- beta.squaredNorm()));
}

double getD(Vector3d &n_los, Vector3d &v) {
    Vector3d beta = v/c;
    return 1./(getG(v)*(1.-beta.dot(n_los)));
}

// (3) from doi:10.1088/0004-637X/695/1/503 (Gracia), (A16) from doi:10.1088/0004-637X/737/1/42 (Porth)
Vector3d get_n_los_prime(Vector3d &n_los, Vector3d &v) {
    double df = getD(n_los, v);
    double gamma = getG(v);
    Vector3d beta = v/c;
    return df*n_los - (df+1.)*(gamma/(gamma+1.))*beta;
}

// (4) in Lyutikov 2003
//Vector3d get_n_los_prime(Vector3d &n_los, Vector3d &v) {
//    double gamma = getG(v);
//    Vector3d beta = v/c;
//    return (n_los + gamma*beta*(gamma*n_los.dot(beta)/(gamma+1)-1))/(gamma*(1-n_los.dot(beta)));
//}

// (4) from doi:10.1088/0004-637X/695/1/503 (Gracia); (4.33) from Daniel, Herbert (1997), "Physik: Elektrodynamik,
// relativistische Physik"
// TODO: Formally, for b = b1 + b2, get_B_prime(b, v) = get_B_prime(b1, v) + get_B_prime(b2, v). Is it OK physically?
Vector3d get_B_prime(Vector3d &b, Vector3d &v) {
    double gamma = getG(v);
    Vector3d result = b/gamma + gamma/((1.+gamma)*(c*c))*v*v.dot(b);
    return result;
}


Vector3d get_B_hat(Vector3d &b_hat_prime, Vector3d &v_beta) {
    double Gamma = 1.0/sqrt(1.0-v_beta.squaredNorm());
    double bv = b_hat_prime.dot(v_beta);
    return (1./sqrt(1.0-bv*bv))*(b_hat_prime - (Gamma/(Gamma+1.0))*bv*v_beta);
}


// from doi:10.1088/0004-637X/737/1/42 (Porth)
//Vector3d get_B_prime(Vector3d &b, Vector3d &v) {
//    double gamma = getG(v);
//    Vector3d result = b/gamma + gamma*v*v.dot(b)/(c*c);
//    return result;
//}

Vector3d q(Vector3d &b_hat, Vector3d &v_beta, Vector3d &n_hat) {
    if (v_beta.norm() == 0){
        return b_hat;
    } else {
        return b_hat + n_hat.cross(v_beta.cross(b_hat));
    }
}

Vector3d e_hat(Vector3d &b_hat, Vector3d &v_beta, Vector3d &n_hat) {
    auto q_ = q(b_hat, v_beta, n_hat);
    auto nq = n_hat.dot(q_);
    return n_hat.cross(q_)/sqrt(q_.squaredNorm() - nq*nq);
}

// This uses adaptive integration as the Astropy version with the same tolerances.
double comoving_transfer_distance2(double z, double H0, double omega_M, double omega_V, double gamma_nu) {
    Ctd ctd(z, H0, omega_M, omega_V, gamma_nu);
    double ctd_state = 0.0;

    using namespace boost::numeric::odeint;
    typedef runge_kutta_dopri5<double> stepper_type;
    auto stepper = stepper_type();
    integrate_adaptive(make_controlled(1.49e-8, 1.49e-8, stepper_type()), ctd, ctd_state, 0.0, z, 1E-09);

    double result = (299792.458/H0) * pow(10., 6.) * ctd_state;
    return result;
}

double pc_to_mas(double z) {
	double d_a = comoving_transfer_distance2(z)/(1.+z);
	double angle_rads = 1./d_a;
	return rad_to_mas*angle_rads;
}

double mas_to_pc(double z) {
	double d_a = comoving_transfer_distance2(z)/(1.+z);
	return mas_to_rad*d_a;
}

std::ostream &
write_2dvector(std::ostream &os, std::vector<std::vector<double>> &v,
							 double scale) {
	for (int i = 0; i < v.size(); ++i) {
		for (int j = 0; j < v[i].size(); ++j) {
			double value = v[i][j]/scale;
			os << value <<" ";
		}
		os<<"\n";
	}
	return os;
}

Ctd::Ctd(double z, double H0, double omega_M, double omega_V, double gamma_nu) : z(z), H0(H0), omega_M(omega_M), omega_V(omega_V), gamma_nu(gamma_nu) {}

void Ctd::operator()(const double &x, double &dxdt, const double t) {
	dxdt = pow((omega_M + (1.+t)*gamma_nu)*(1.+t)*(1.+t)*(1.+t)+omega_V, -0.5);
	//dxdt = pow(omega_M*(1.+t*t*t)+omega_V, -0.5);
};


void read_from_txt(const std::string& fntxt, std::vector< std::vector<double> >& properties) {
    std::ifstream infile(fntxt);
    if(!infile.good()){
        throw AbsentDataFile(fntxt);
    }
    std::vector<double> row(3);
    while (infile >> row[0] >> row[1] >> row[2])
    {
        properties.push_back(row);
    }
}
