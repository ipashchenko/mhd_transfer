#include <utils.h>
#include "NField.h"
#include "MyExceptions.h"
//#include "symphony.h"


NField::NField (bool in_plasma_frame, Geometry* geometry) :
        in_plasma_frame_(in_plasma_frame) {
    geometry_ = geometry;
}


double NField::nf(const Vector3d &point, double psi) const {
    double x, y, r_point, r_border, result;

    if(geometry_) {
        // Find radius at given point
        r_border = geometry_->radius_at_given_distance(point);
        x = point[0];
        y = point[1];
        r_point = sqrt(x*x + y*y);
        if (r_point > r_border) {
            return 0.0;
        }
        else {
            result = _nf(point, psi);
        }
    }
    else {
        result = _nf(point, psi);
    }

    return result;
}

double NField::nf_plasma_frame(const Vector3d &point, double psi, double &gamma) const {
    double n = nf(point, psi);
    if (in_plasma_frame_) {
        return n;
    } else {
        return n/gamma;
    }
}


double NField::k_F_c(Vector3d &b, Vector3d &n_los, double nu, double n) const {
    double cos_theta = b.dot(n_los)/b.norm();
    return 2.*pi*nu_p(n)*nu_p(n)*nu_b_value(b)*cos_theta/(c*nu*nu);
}

double NField::k_C_c(Vector3d &b, Vector3d &n_los, double nu, double n) const {
    auto nu_b_calc = nu_b(b, n_los);
    return -pi*nu_p(n)*nu_p(n)*nu_b_calc*nu_b_calc/(c*nu*nu*nu);
}


PowerLawNField:: PowerLawNField(bool in_plasma_frame, double s, double gamma_min, Geometry* geometry,
                               std::string plasma, bool changing_s, double ds) :
        NField(in_plasma_frame, geometry), s_(s), gamma_min_(gamma_min), plasma_(plasma), changing_s_(changing_s), ds_(ds) {
    if(plasma_ == "normal"){
        std::cout << "Normal plasma content" << std::endl;
    }
    else if(plasma_ == "pairs"){
        std::cout << "Pairs plasma content" << std::endl;
    } else {
        throw BadPlasmaContentParameters();
    }
    factor_ki_ = (pow(3., (s_+1.)/2.)/4.)*tgamma(s_/4.+11./6.)*tgamma(s_/4.+1./6.);
    factor_ki_rnd_ = sqrt(pi/4.)*tgamma((6.+s_)/4.)/tgamma((8.+s_)/4.);
    factor_kv_ = pow(3., s_/2.)*(s_+3.)*(s_+2.)/(4.*(s_+1.))*tgamma(s_/4.+11./12.)*tgamma(s_/4.+7./12.);
    factor_kf_ = (s_+2.)*log(gamma_min_)/((s_+1.)*pow(gamma_min_, s_+1.));
    factor_etai_ = pow(3., s_/2.)/(2.*(s_+1))*tgamma(s_/4.+19./12.)*tgamma(s_/4.-1./12.);
    factor_etai_rnd_ = sqrt(pi/4.)*tgamma((5.+s_)/4.)/tgamma((7.+s_)/4.);
    factor_etav_ = pow(3., (s_-1.)/2.)*(s_+2.)/(2.*s_)*tgamma(s_/4.+2./3.)*tgamma(s_/4.+1./3.);
    if(plasma_ == "pairs"){
        factor_kf_ = 0.0;
        factor_kv_ = 0.0;
        factor_etav_ = 0.0;
    }
}

double PowerLawNField::get_s(const Vector3d &point) const {
    double r = point.norm();
    return s_ + ds_*r/pc;
}

double PowerLawNField::get_s(Vector3d &b, Vector3d &n_los) const {
    double s_along = 2;
    double s_across = 2.5;
    double cos_theta = abs(b.dot(n_los)/b.norm());
//    std::cout << "s = " << s_along*cos_theta + s_across*(1-cos_theta) << "\n";
    return s_along*cos_theta + s_across*(1-cos_theta);
}

double PowerLawNField::k_i(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const {
    double factor_ki;
    double s;
    if (changing_s_) {
        s = get_s(b, n_los);
        factor_ki = (pow(3., (s+1.)/2.)/4.)*tgamma(s/4.+11./6.)*tgamma(s/4.+1./6.);
    } else {
        s = s_;
        factor_ki = factor_ki_;
    }
    double nu_min_ = nu_min(gamma_min_, b, n_los);
    if(nu > nu_min_) {
        return k_0(b, n_los, nu, n_nt, s, gamma_min_) * pow(nu_b(b, n_los)/nu, s/2.) * factor_ki;
//        return k_0(b, n_los, nu, n*(s-1)*pow(gamma_min_, s-1)) * pow(nu_b(b, n_los)/nu, s/2.) * factor_ki;
    } else {
        return k_0(b, n_los, nu_min_, n_nt, s, gamma_min_) * pow(nu_b(b, n_los)/nu_min_, s/2.) * factor_ki * pow(nu/nu_min_, -5./3.);
//        return k_0(b, n_los, nu_min_, n*(s-1)*pow(gamma_min_, s-1)) * pow(nu_b(b, n_los)/nu_min_, s/2.) * factor_ki * pow(nu/nu_min_, -5./3.);
    }
}

double PowerLawNField::k_i(double b, Vector3d &n_los, double nu, double n_nt) const {
    //double factor = (pow(3., (s_+1.)/2.)/4.)*tgamma(s_/4.+11./6.)*tgamma(s_/4.+1./6.);
    //double rnd_factor = sqrt(pi/4.)*tgamma((6.+s_)/4.)/tgamma((8.+s_)/4.);
    double factor = factor_ki_*factor_ki_rnd_;
    return k_0(b, n_los, nu, n_nt, s_, gamma_min_) * pow(nu_b(b)/nu, s_/2.) * factor;
}

double PowerLawNField::k_q(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const {
    double s;
    if (changing_s_) {
        s = get_s(b, n_los);
    } else {
        s = s_;
    }
    return (s+2.)/(s+10./3)*k_i(b, n_los, nu, n_nt);
}

double PowerLawNField::k_u(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const {
    return 0;
}

double PowerLawNField::k_v(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const {

    double factor_kv;
    double s;
    if (changing_s_) {
        s = get_s(b, n_los);
        if(plasma_ == "pairs"){
            factor_kv = 0;
        } else {
            factor_kv = pow(3., s/2.)*(s+3.)*(s+2.)/(4.*(s+1.))*tgamma(s/4.+11./12.)*tgamma(s/4.+7./12.);
        }
    } else {
        s = s_;
        factor_kv = factor_kv_;
    }

    double cos_theta = b.dot(n_los)/b.norm();
    //double factor = pow(3., s_/2.)*(s_+3.)*(s_+2.)/(4.*(s_+1.))*tgamma(s_/4.+11./12.)*tgamma(s_/4.+7./12.);
    double nu_min_ = nu_min(gamma_min_, b, n_los);
    if(nu > nu_min_) {
        return -k_0_value(b, nu, n_nt, s, gamma_min_)*cos_theta*pow(nu_b(b, n_los)/nu, (s+1.)/2.) * factor_kv;
    } else {
        return -k_0_value(b, nu_min_, n_nt, s, gamma_min_)*cos_theta*pow(nu_b(b, n_los)/nu_min_, (s+1.)/2.) * factor_kv * pow(nu/nu_min_, -5./3.);
    }
}

double PowerLawNField::k_F(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const {
    double factor_kf;
    double s;
    if (changing_s_) {
        s = get_s(b, n_los);
        if(plasma_ == "pairs") {
            factor_kf = 0.0;
        } else {
            factor_kf = (s+2.)*log(gamma_min_)/((s+1.)*pow(gamma_min_, s+1.));
        }
    } else {
        factor_kf = factor_kf_;
    }
    //return (s_+2.)*log(gamma_min_)/((s_+1.)*pow(gamma_min_, s_+1.))*k_F_c(b, n_los, nu, n);
    return factor_kf * k_F_c(b, n_los, nu, n_nt);
}

double PowerLawNField::k_C(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const {
    double s;
    if (changing_s_) {
        s = get_s(b, n_los);
    } else {
        s = s_;
    }
    return (2./(s-2.))*(pow(gamma_min_, 2.-s) - pow(nu_b(b, n_los)/nu, (s-2.)/2.))*k_C_c(b, n_los, nu, n_nt);
}

double PowerLawNField::h_Q(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const {
    return 0;
}

double PowerLawNField::eta_i(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const {
    //double factor = pow(3., s_/2.)/(2.*(s_+1))*tgamma(s_/4.+19./12.)*tgamma(s_/4.-1./12.);

    double factor_etai;
    double s;
    if (changing_s_) {
        s = get_s(b, n_los);
        factor_etai = pow(3., s/2.)/(2.*(s+1))*tgamma(s/4.+19./12.)*tgamma(s/4.-1./12.);
    } else {
        s = s_;
        factor_etai = factor_etai_;
    }

    double nu_min_ = nu_min(gamma_min_, b, n_los);
    if(nu > nu_min_) {
        return eta_0(b, n_los, n_nt, s, gamma_min_) * pow(nu_b(b, n_los)/nu, (s-1.)/2.) * factor_etai;
//        return eta_0(b, n_los, n*(s-1)*pow(gamma_min_, s-1)) * pow(nu_b(b, n_los)/nu, (s-1.)/2.) * factor_etai;
    } else {
        return eta_0(b, n_los, n_nt, s, gamma_min_) * pow(nu_b(b, n_los)/nu_min_, (s-1.)/2.) * factor_etai * pow(nu/nu_min_, 1./3.);
//        return eta_0(b, n_los, n*(s-1)*pow(gamma_min_, s-1)) * pow(nu_b(b, n_los)/nu_min_, (s-1.)/2.) * factor_etai * pow(nu/nu_min_, 1./3.);
    }
}

double PowerLawNField::eta_i(double b, Vector3d &n_los, double nu, double n_nt) const {
    //double factor = pow(3., s_/2.)/(2.*(s_+1))*tgamma(s_/4.+19./12.)*tgamma(s_/4.-1./12.);
    //double rnd_factor = sqrt(pi/4.)*tgamma((5.+s_)/4.)/tgamma((7.+s_)/4.);
    double factor = factor_etai_ * factor_etai_rnd_;
    return eta_0(b, n_los, n_nt, s_, gamma_min_) * pow(nu_b(b)/nu, (s_-1.)/2.) * factor;
}

double PowerLawNField::eta_q(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const {
    double s;
    if (changing_s_) {
        s = get_s(b, n_los);
    } else {
        s = s_;
    }
    return (s+1.0)/(s+7./3.)*eta_i(b, n_los, nu, n_nt);
}

double PowerLawNField::eta_u(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const {
    return 0;
}

double PowerLawNField::eta_v(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const {

    double factor_etav;
    double s;
    if (changing_s_) {
        s = get_s(b, n_los);
        if(plasma_ == "pairs") {
            factor_etav = 0.0;
        } else {
            factor_etav = pow(3., (s-1.)/2.)*(s+2.)/(2.*s)*tgamma(s/4.+2./3.)*tgamma(s/4.+1./3.);
        }
    } else {
        s = s_;
        factor_etav = factor_etav_;
    }

    double cos_theta = b.dot(n_los)/b.norm();
    //double factor = pow(3., (s_-1.)/2.)*(s_+2.)/(2.*s_)*tgamma(s_/4.+2./3.)*tgamma(s_/4.+1./3.);

    double nu_min_ = nu_min(gamma_min_, b, n_los);
    if(nu > nu_min_) {
        return -eta_0_value(b, n_nt, s, gamma_min_)*cos_theta*pow(nu_b(b, n_los)/nu, s/2.) * factor_etav;
    } else {
        return -eta_0_value(b, n_nt, s, gamma_min_)*cos_theta*pow(nu_b(b, n_los)/nu_min_, s/2.) * factor_etav * pow(nu/nu_min_, 1./3.);
    }
}


ThermalNField::ThermalNField(Geometry* geometry) : NField(true, geometry) {}

double ThermalNField::k_i(Vector3d &b, Vector3d &n_los, double nu, double n) const {
    return 0.0;
}

double ThermalNField::k_i(double b, Vector3d &n_los, double nu, double n) const {
    return 0.0;
}

double ThermalNField::k_q(Vector3d &b, Vector3d &n_los, double nu, double n) const {
    return 0.0;
}

double ThermalNField::k_u(Vector3d &b, Vector3d &n_los, double nu, double n) const {
    return 0.0;
}

double ThermalNField::k_v(Vector3d &b, Vector3d &n_los, double nu, double n) const {
    return 0.0;
}

double ThermalNField::k_F(Vector3d &b, Vector3d &n_los, double nu, double n) const {
    return k_F_c(b, n_los, nu, n);
}

double ThermalNField::k_C(Vector3d &b, Vector3d &n_los, double nu, double n) const {
    return k_C_c(b, n_los, nu, n);
}

double ThermalNField::h_Q(Vector3d &b, Vector3d &n_los, double nu, double n) const {
    return 0.0;
}

double ThermalNField::eta_i(Vector3d &b, Vector3d &n_los, double nu, double n) const {
    return 0.0;
}

double ThermalNField::eta_i(double b, Vector3d &n_los, double nu, double n) const {
    return 0.0;
}

double ThermalNField::eta_q(Vector3d &b, Vector3d &n_los, double nu, double n) const {
    return 0.0;
}

double ThermalNField::eta_u(Vector3d &b, Vector3d &n_los, double nu, double n) const {
    return 0.0;
}

double ThermalNField::eta_v(Vector3d &b, Vector3d &n_los, double nu, double n) const {
    return 0.0;
}


SimulationNField::SimulationNField(Delaunay_triangulation *tr, bool in_plasma_frame, double s, double gamma_min,
                                   bool changing_s, double scale_factor) :
        PowerLawNField(in_plasma_frame, s, gamma_min, nullptr, "pairs", changing_s, 0.0),
        interp_(tr),
        scale_factor_(scale_factor) {}

double SimulationNField::_nf(const Vector3d &point, double psi) const {
    double z = point[2];
    return scale_factor_*interp_.interpolated_value({psi, z/pc});
};


ConstrainedSimulationNField::ConstrainedSimulationNField(Delaunay_triangulation *tr_nt,
                                                         Delaunay_triangulation *tr_cold,
                                                         bool in_plasma_frame, double s, double gamma_min,
                                                         bool changing_s, double scale_factor_nt,
                                                         double scale_factor_cold) :
        PowerLawNField(in_plasma_frame, s, gamma_min, nullptr, "pairs", changing_s, 0.0),
        interp_nt_(tr_nt),
        interp_cold_(tr_cold),
        scale_factor_nt_(scale_factor_nt),
        scale_factor_cold_(scale_factor_cold) {}

double ConstrainedSimulationNField::_nf(const Vector3d &point, double psi) const {
    double z = point[2];

    double n_cold = scale_factor_cold_*interp_cold_.interpolated_value({psi, z/pc});
    double n_nt = scale_factor_nt_*interp_nt_.interpolated_value({psi, z/pc});
    if(n_nt < n_cold) {
        return n_nt;
    }else {
        return n_cold;
    }
};


ConstrainedSigmaSimulationNField::ConstrainedSigmaSimulationNField(Delaunay_triangulation *tr_nt,
                                                                   Delaunay_triangulation *tr_cold,
                                                                   Delaunay_triangulation *tr_sigma,
                                                                   bool in_plasma_frame, double s, double gamma_min,
                                                                   bool changing_s, double scale_factor_nt,
                                                                   double scale_factor_cold) :
        PowerLawNField(in_plasma_frame, s, gamma_min, nullptr, "pairs", changing_s, 0.0),
        interp_nt_(tr_nt),
        interp_cold_(tr_cold),
        interp_sigma_(tr_sigma),
        scale_factor_nt_(scale_factor_nt),
        scale_factor_cold_(scale_factor_cold) {}

double ConstrainedSigmaSimulationNField::_nf(const Vector3d &point, double psi) const {
    double z = point[2];

    double sigma = interp_sigma_.interpolated_value({psi, z/pc});
    double sigma_suppression_factor = 1./(1. + exp(sigma));

    double n_cold = scale_factor_cold_*interp_cold_.interpolated_value({psi, z/pc});
    double n_nt = sigma_suppression_factor*scale_factor_nt_*interp_nt_.interpolated_value({psi, z/pc});
    if(n_nt < n_cold) {
        return n_nt;
    }else {
        return n_cold;
    }
};


ConstrainedBetaSimulationNField::ConstrainedBetaSimulationNField(Delaunay_triangulation *tr_cold,
                                                                 Delaunay_triangulation *tr_Bsq,
                                                                 Delaunay_triangulation *tr_jsq,
                                                                 Delaunay_triangulation *tr_sigma,
                                                                 std::string particles_heating_model,
                                                                 std::string constrain_type,
                                                                 bool in_plasma_frame, double s, double gamma_min,
                                                                 bool changing_s,
                                                                 double scale_factor_nt,
                                                                 double max_frac_cold,
                                                                 double scale_factor_border,
                                                                 double psi_mean,
                                                                 double psi_width) :
        PowerLawNField(in_plasma_frame, s, gamma_min, nullptr, "pairs", changing_s, 0.0),
        interp_cold_(tr_cold),
        interp_Bsq_(tr_Bsq),
        interp_jsq_(tr_jsq),
        interp_sigma_(tr_sigma),
        particles_heating_model_(std::move(particles_heating_model)),
        constrain_type_(std::move(constrain_type)),
        scale_factor_nt_(scale_factor_nt),
        scale_factor_border_(scale_factor_border),
        max_frac_cold_(max_frac_cold),
        psi_mean_(psi_mean),
        psi_width_(psi_width) {}

double ConstrainedBetaSimulationNField::_nf(const Vector3d &point, double psi) const {
    double x = point[0];
    double y = point[1];
    double r_p = hypot(x, y)/pc;
    double z = point[2];
    double n_nt;


    // Relative amplitude of NT particles density enhancement
//    double amp = -(psi_mean_ - psi)*(psi_mean_ - psi)/(2.0*psi_width_*psi_width_);
//    double amp = -abs(psi_mean_ - psi)/psi_width_;
//    amp = exp(amp);

//    double amp;
//    if(abs(psi - psi_mean_) < psi_width_){
//        amp = 1.0;
//    }else {
//        amp = 0.0;
//    }

    double amp = generalized1_gaussian1d(psi, psi_mean_, psi_width_, 2);

    // Number of cold particles in plasma frame
    double n_cold = interp_cold_.interpolated_value({psi, z/pc});
    // Maximum number of cold particles that can be heated
    double n_cold_max = max_frac_cold_*n_cold;

    // Squared magnetic field in plasma frame
    double Bsq;
    if(constrain_type_ == "beta" || particles_heating_model_ == "bsq") {
        Bsq = interp_Bsq_.interpolated_value({psi, z/pc});
    }


    // NT particles heating can be reduced in some regions
    double suppression_factor;
    if(constrain_type_ == "beta") {
        // Local inverse sigma
        double inv_beta = 8.0 * M_PI * n_cold * m_e * c * c / Bsq;
        // Suppression factor of Broderick+2010
        suppression_factor = 1./(1. + exp(inv_beta));
    }
    else if(constrain_type_ == "sigma") {
        double sigma = interp_sigma_.interpolated_value({psi, z/pc});
        // Lena's suggestion
        suppression_factor = 1./(1. + exp(1./sigma));
    }
    else if(constrain_type_ == "none") {
        // No modifications
        suppression_factor = 1.0;
    }
    else {
        throw NotImplmentedHeatingSuppression(constrain_type_);
    }

    double n_nt_border;

    if(particles_heating_model_ == "n"){
        // Number of emitting particles ~ number of cold particles
        // In this case ``scale_factor_nt_`` - efficiency of cold particles heating - must be in [0, 1]
//        n_nt = scale_factor_nt_*n_cold*suppression_factor;
        n_nt = n_cold*suppression_factor;
    }
    else if(particles_heating_model_ == "bsq") {
        // Number of emitting particles ~ B_plasma^2
        // In this case ``scale_factor_nt_`` must be in [0, 1]. It is efficiency with which magnetic energy is converted
        // into that of non-thermal particles.
        // From the relation:
        // u_e = n_nt * mc^2 * (s-1)/(s-2) * gamma_min = B^2/(8pi) * sigma_suppression_factor
        // the scale coefficient from B^2 to n_nt is:
        double n_nt_Bsq = suppression_factor*(s_ - 2)/(s_ - 1)/(8*M_PI*m_e*c*c*gamma_min_);
//        n_nt = scale_factor_nt_*n_nt_Bsq*Bsq;
        n_nt = n_nt_Bsq*Bsq;
    }
    else if(particles_heating_model_ == "jsq") {
        // Number of emitting particles ~ j_plasma^2. ``scale_factor_nt_`` is arbitrary positive number.
//        n_nt = scale_factor_nt_*interp_jsq_.interpolated_value({psi, z/pc})*suppression_factor;
        n_nt = interp_jsq_.interpolated_value({psi, z/pc})*suppression_factor;
    }
    // FIXME: This is rather unrealistic because NT particles density doesn't change with distance along the jet.
    else if(particles_heating_model_ == "byhand") {
        n_nt = 1.0;
    }
    else {
        throw NotImplmentedParticlesHeating(particles_heating_model_);
    }

    // Artificial humps:
    n_nt_border = scale_factor_border_ * amp * n_nt;
    // Model heating:
    // From model heating. Here ``scale_factor_nt_`` must be [0, 1] for heating models ``n``, ``bsq``; any positive
    // number for ``jsq``, ``byhand``.
    n_nt = scale_factor_nt_ * n_nt;

    // Number of non-thermal particles can't be larger than some fraction (``max_frac_cold_``) of all cold particles
    if(n_nt + n_nt_border > n_cold_max) {
//        std::cout << "Number of NT particles can't be larger than number of cold particles!" << std::endl;
        return n_cold_max;
    }else {
        return n_nt + n_nt_border;
    }
};


ByHandSimulationNField::ByHandSimulationNField(Delaunay_triangulation *tr_cold,
                                               bool in_plasma_frame, double s, double gamma_min,
                                               bool changing_s,
                                               double max_frac_cold,
                                               double scale_factor_border,
                                               double psi_mean,
                                               double psi_width,
                                               double psi_width_axis,
                                               double scale_factor_axis,
                                               double constant_floor_scale,
                                               double n) :
        PowerLawNField(in_plasma_frame, s, gamma_min, nullptr, "pairs", changing_s, 0.0),
        interp_cold_(tr_cold),
        scale_factor_border_(scale_factor_border),
        scale_factor_axis_(scale_factor_axis),
        max_frac_cold_(max_frac_cold),
        psi_mean_(psi_mean),
        psi_width_(psi_width),
        psi_width_axis_(psi_width_axis),
        constant_floor_scale_(constant_floor_scale),
        n_(n) {}

double ByHandSimulationNField::_nf(const Vector3d &point, double psi) const {
    double z = point[2];

    // Outer ridge
    // Step profile
//    double amp_border;
//    if(abs(psi - psi_mean_) < psi_width_){
//        amp_border = 1.0;
//    }else {
//        amp_border = 0.0;
//    }
    // Generalized Gaussian profile with maximum equal to 1.0 at ``psi_mean_``
    double amp_border = generalized1_gaussian1d(psi, psi_mean_, psi_width_, 2);
    double amp_axis = generalized1_gaussian1d(psi, 0.0, psi_width_axis_, 2);

    // Number of cold particles in plasma frame at current given point
    double n_cold = interp_cold_.interpolated_value({psi, z/pc});
    double n_nt_border = scale_factor_border_ * amp_border * pow(z / pc, -n_);
    double n_nt_axis = scale_factor_axis_ * amp_axis * pow(z/pc, -n_);
    double n_nt_const = constant_floor_scale_ * pow(z/pc, -n_);

    // Maximum number of cold particles that can be heated
    double n_cold_max = max_frac_cold_*n_cold;

    // Number of non-thermal particles can't be larger than some fraction (``max_frac_cold_``) of all cold particles
    if(n_nt_border + n_nt_axis + n_nt_const > n_cold_max) {
        std::cout << "Number of NT particles can't be larger than number of cold particles!" << std::endl;
        return n_cold_max;
    }else {
        return n_nt_border + n_nt_axis + n_nt_const;
    }
};