#ifndef MHD_TRANSFER_UTILS_H
#define MHD_TRANSFER_UTILS_H

#include <cmath>
#include <vector>
#include <Eigen/Eigen>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/odeint.hpp>
#include "Geometry.h"

using Eigen::Vector3d;


const double mas_to_rad = 4.84813681109536e-09;
const double rad_to_mas = 206264806.24709633;
// Parsec [cm]
const double pc = 3.0856775814671913e+18;
// Mass of electron [g]
const double m_e = 9.10938356e-28;
// Mass of proton [g]
const double m_p = 1.672621898e-24;
// Charge of electron [C]
 const double q_e = 4.80320467299766e-10;
//const double q_e = 1.6*1E-19;
// Charge of proton [C]
 const double q_p = q_e;
//const double q_p = 1.6*1E-19;
// Speed of light [cm / s]
const double c = 29979245800.0;
// Gravitational constant
const double G = 6.67430e-08;

// Jy in cgc
const double to_jy = 1E+23;
// pi
const double pi = boost::math::constants::pi<double>();

// Minimal value of B to calculate transfer coefficients
const double eps_B = 1E-10;


double nu_p(double n);

double nu_b(double b);

double nu_b(Vector3d &b, Vector3d &n_los);

double nu_b_value(Vector3d &b);

double sin_theta(Vector3d &b, Vector3d &n_los);

// Frequency that corresponds to lower energy cut-off gamma_min in emitting electrons energy spectrum
double nu_min(double gamma_min, Vector3d &b, Vector3d &n_los);

double nu_min(double gamma_min, double b);

// For random B-field
double k_0(double b, Vector3d &n_los, double nu, double n);

// For vector B-field
double k_0(Vector3d &b, Vector3d &n_los, double nu, double n);

double k_0_value(Vector3d &b, double nu, double n);

// For random B-field
double k_i(double b, Vector3d &n_los, double nu, double n, double s=2.5);

// For vector B-field
// Absorption coefficient for given vector of magnetic field ``b``, unit LOS
// vector ``n_los`` and others measured in emission frame
double k_i(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double eta_0(double b, Vector3d &n_los, double n);

// For vector B-field
double eta_0(Vector3d &b, Vector3d &n_los, double n);

double eta_0_value(Vector3d &b, double n);

// For random B-field
double eta_i(double b, Vector3d &n_los, double nu, double n, double s=2.5);

// For vector B-field
// Emission coefficient for given vector of magnetic field ``b``, unit LOS
// vector ``n_los`` and others measured in emission frame
double eta_i(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For vector B-field
double k_q(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double k_q(double b, Vector3d &n_los, double nu, double n, double s=2.5);

// For vector B-field
double k_u(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double k_u(double b, Vector3d &n_los, double nu, double n, double s=2.5);

// For vector B-field
double k_v(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double k_v(double b, Vector3d &n_los, double nu, double n, double s=2.5);

// For vector B-field
double k_F_c(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double k_F_c(double b, Vector3d &n_los, double nu, double n, double s=2.5);

// For vector B-field
double k_C_c(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double k_C_c(double b, Vector3d &n_los, double nu, double n, double s=2.5);

// For vector B-field
// TODO: For e+ ``k_F_r`` will have different sign, but ``k_C_r`` - the same
// sign
double k_F_r(Vector3d &b, Vector3d &n_los, double nu, double n,
             double gamma_min=100., double s=2.5);

// For random B-field
double k_F_r(double b, Vector3d &n_los, double nu, double n,
             double gamma_min=100., double s=2.5);

// For vector B-field
double k_C_r(Vector3d &b, Vector3d &n_los, double nu, double n,
             double gamma_min=100., double s=2.5);

// For random B-field
double k_C_r(double b, Vector3d &n_los, double nu, double n,
             double gamma_min=100., double s=2.5);

// For vector B-field
double h_Q(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double h_Q(double b, Vector3d &n_los, double nu, double n, double s=2.5);

// For vector B-field
double eta_q(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double eta_q(double b, Vector3d &n_los, double nu, double n, double s=2.5);

double eta_u(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double eta_u(double b, Vector3d &n_los, double nu, double n, double s=2.5);

double eta_v(Vector3d &b, Vector3d &n_los, double nu, double n, double s=2.5);

// For random B-field
double eta_v(double b, Vector3d &n_los, double nu, double n, double s=2.5);

// Lorentz factor for the velocity ``v``
double getG(Vector3d &v);

// Doppler factor for LOS unit vector ``n_los`` in lab frame and bulk velocity
// ``v`` relative to the lab frame.
double getD(Vector3d &n_los, Vector3d &v);

// LOS unit vector in the comoving frame. This is relativistic aberration.
Vector3d get_n_los_prime(Vector3d &n_los, Vector3d &v);

// Vector of magnetic field in comoving (emission) frame in terms of the
// magnetic field ``b`` and velocity ``v`` in the observer frame.
Vector3d get_B_prime(Vector3d &b, Vector3d &v);

// Unit vector of magnetic field in laboratory (observer) frame in terms of the unit vector of the magnetic field
// ``b_hat_prime`` in plasma (comoving) and velocity ``v_beta`` (in c units) in the observer frame.
Vector3d get_B_hat(Vector3d &b_hat_prime, Vector3d &v_beta);

// Following two are (C5) from Lyutikov 2003. Needed to express the unit vector of electric field of the wave in the
// observer frame. b_hat - unit vector of B in lab frame, v_beta - vector of speed (in c units) in lab frame, n_hat -
// direction of the emitted photons in lab frame.
Vector3d q(Vector3d &b_hat, Vector3d &v_beta, Vector3d &n_hat);

Vector3d e_hat(Vector3d &b_hat, Vector3d &v_beta, Vector3d &n_hat);


class Ctd {
public:
	//From astropy.cosmology.WMAP9
	Ctd(double z, double H0=69.32, double omega_M=0.2865, double omega_V=0.7134130719051658,
			double gamma_nu=8.69280948342326e-05);
	void operator() (const double &x, double &dxdt, const double t);
	double z;
	double H0;
	double omega_M;
	double omega_V;
	double gamma_nu;
};

double comoving_transfer_distance2(double z, double H0=69.32, double omega_M=0.2865, double omega_V=0.7134130719051658, double gamma_nu=8.69280948342326e-05);

// Return scale factor that converts from parsecs to milliarcseconds
double pc_to_mas(double z);

// Return scale factor that converts from milliarcseconds to parsecs.
double mas_to_pc(double z);

std::ostream& write_2dvector(std::ostream& os, std::vector<std::vector<double>>& v, double scale=1.0);

void read_from_txt(std::string fn, std::vector< std::vector<double> >& properties);

#endif //MHD_TRANSER_UTILS_H
