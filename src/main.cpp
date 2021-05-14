#include <string>
#include <iostream>
#include <fstream>
#include "ImagePlane.h"
#include "Observation.h"
#include "BField.h"
#include "VField.h"
#include "Jet.h"
#include "utils.h"
#include <cmath>
#include "NField.h"
#include "Pixel.h"
#include <ctime>
#include <chrono>
#include "linspace.h"
#include "MyExceptions.h"

using Eigen::Vector3d;
using Eigen::Matrix3Xd;
using std::vector;
using std::pair;
using namespace boost::numeric::odeint;
namespace ph = std::placeholders;

// For running w/o simulations
typedef std::chrono::high_resolution_clock Clock;


void check_interpolations(const std::string& mhd_run_name) {

    Delaunay_triangulation tr_B_p, tr_B_phi, tr_N, tr_Gamma, tr_beta_phi, tr_jsq, tr_Psi;
    create_triangulation(mhd_run_name + "_B_p_field.txt", &tr_B_p);
    create_triangulation(mhd_run_name + "_B_phi_field.txt", &tr_B_phi);
    create_triangulation(mhd_run_name + "_Psi_field.txt", &tr_Psi);
    create_triangulation(mhd_run_name + "_n_plasma_field.txt", &tr_N);
    create_triangulation(mhd_run_name + "_Gamma_field.txt", &tr_Gamma);
    create_triangulation(mhd_run_name + "_beta_phi_field.txt", &tr_beta_phi);
    create_triangulation(mhd_run_name + "_jsq_plasma_field.txt", &tr_jsq);


    SimulationInterpolater interp_B_p(&tr_B_p);
    SimulationInterpolater interp_B_phi(&tr_B_phi);
    SimulationInterpolater interp_N(&tr_N);
    SimulationInterpolater interp_Gamma(&tr_Gamma, 1.0);
    SimulationInterpolater interp_beta_phi(&tr_beta_phi);
    SimulationInterpolater interp_jsq(&tr_jsq);
    SimulationInterpolater interp_Psi(&tr_Psi, 1.0);

    // Points where to find values
    size_t n_along = 1500;
    size_t n_across = 200;
    auto z_pc = linspace(0.0, 1.5, 1501);
    auto r_pc = linspace(0.0, 0.2, 201);

    std::vector<std::vector<double>> B_p, B_phi, N, Gamma, beta_phi, jsq, Psi;
    B_p.resize(n_along);
    B_phi.resize(n_along);
    N.resize(n_along);
    Gamma.resize(n_along);
    beta_phi.resize(n_along);
    jsq.resize(n_along);
    Psi.resize(n_along);
    for(size_t i=0; i < n_along; i++) {
        B_p[i].resize(n_across);
        B_phi[i].resize(n_across);
        N[i].resize(n_across);
        Gamma[i].resize(n_across);
        beta_phi[i].resize(n_across);
        jsq[i].resize(n_across);
        Psi[i].resize(n_across);
        for(size_t j=0; j < n_across; j++) {
            Vector2d pos = {r_pc[j], z_pc[i]};
            B_p[i][j] = interp_B_p.interpolated_value(pos);
            B_phi[i][j] = interp_B_phi.interpolated_value(pos);
            N[i][j] = interp_N.interpolated_value(pos);
            Gamma[i][j] = interp_Gamma.interpolated_value(pos);
            beta_phi[i][j] = interp_beta_phi.interpolated_value(pos);
            jsq[i][j] = interp_jsq.interpolated_value(pos);
            Psi[i][j] = interp_Psi.interpolated_value(pos);
        }
    }

    std::fstream fs;

    fs.open(mhd_run_name + "_B_p_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, B_p);
        fs.close();
    }

    fs.open(mhd_run_name + "_B_phi_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, B_phi);
        fs.close();
    }

    fs.open(mhd_run_name + "_N_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, N);
        fs.close();
    }

    fs.open(mhd_run_name + "_Gamma_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, Gamma);
        fs.close();
    }

    fs.open(mhd_run_name + "_beta_phi_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, beta_phi);
        fs.close();
    }

    fs.open(mhd_run_name + "_jsq_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, jsq);
        fs.close();
    }

    fs.open(mhd_run_name + "_Psi_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, Psi);
        fs.close();
    }
}


void check_psi_interpolations_Lena(const std::string& mhd_run_name) {

    Delaunay_triangulation tr_B_p, tr_B_phi, tr_N, tr_Gamma, tr_beta_phi, tr_jsq, tr_jsq_phi, tr_sigma, tr_theta, tr_Psi, tr_rPsi;
    // Triangulations for (z, Psi, value)
    create_triangulation(mhd_run_name + "_B_p_field_psi.txt", &tr_B_p);
    create_triangulation(mhd_run_name + "_B_phi_field_psi.txt", &tr_B_phi);
    create_triangulation(mhd_run_name + "_Psi_field_psi.txt", &tr_Psi);
    create_triangulation(mhd_run_name + "_n_plasma_field_psi.txt", &tr_N);
    create_triangulation(mhd_run_name + "_Gamma_field_psi.txt", &tr_Gamma);
    create_triangulation(mhd_run_name + "_beta_phi_field_psi.txt", &tr_beta_phi);
    create_triangulation(mhd_run_name + "_sigma_field_psi.txt", &tr_sigma);
    create_triangulation(mhd_run_name + "_theta_field_psi.txt", &tr_theta);
    create_triangulation(mhd_run_name + "_jsq_phi_plasma_field_psi.txt", &tr_jsq_phi);
    create_triangulation(mhd_run_name + "_jsq_plasma_field_psi.txt", &tr_jsq);
    // Triangulation for (z, r, Psi)
    create_triangulation(mhd_run_name + "_Psi_field.txt", &tr_rPsi);

    // Interpolates in (z, Psi) coordinates
    SimulationInterpolater interp_B_p(&tr_B_p);
    SimulationInterpolater interp_B_phi(&tr_B_phi);
    SimulationInterpolater interp_N(&tr_N);
    SimulationInterpolater interp_Gamma(&tr_Gamma, 1.0);
    SimulationInterpolater interp_sigma(&tr_sigma);
    SimulationInterpolater interp_theta(&tr_theta);
    SimulationInterpolater interp_beta_phi(&tr_beta_phi);
    SimulationInterpolater interp_jsq(&tr_jsq);
    SimulationInterpolater interp_jsq_phi(&tr_jsq_phi);
    SimulationInterpolater interp_Psi(&tr_Psi, 1.0);

    // Interpolate Psi(z, r): given (z, r) returns (Psi)
    SimulationInterpolater interp_rPsi(&tr_rPsi, 1.0);
    // Points where to find values
    size_t n_along = 5000;
    size_t n_across = 500;
    auto z_pc = linspace(0.0, 5.0, 5001);
    auto r_pc = linspace(0.0, 0.5, 501);
    auto psi = linspace(0.0, 1.0, 501);

    std::vector<std::vector<double>> B_p, B_phi, N, Gamma, beta_phi, jsq, jsq_phi, sigma, theta, Psi, angle;
    B_p.resize(n_along);
    B_phi.resize(n_along);
    N.resize(n_along);
    Gamma.resize(n_along);
    beta_phi.resize(n_along);
    jsq.resize(n_along);
    jsq_phi.resize(n_along);
    sigma.resize(n_along);
    theta.resize(n_along);
    Psi.resize(n_along);
    angle.resize(n_along);
    for(size_t i=0; i < n_along; i++) {
        B_p[i].resize(n_across);
        B_phi[i].resize(n_across);
        N[i].resize(n_across);
        Gamma[i].resize(n_across);
        beta_phi[i].resize(n_across);
        jsq[i].resize(n_across);
        jsq_phi[i].resize(n_across);
        sigma[i].resize(n_across);
        theta[i].resize(n_across);
        Psi[i].resize(n_across);
        angle[i].resize(n_across);
        for(size_t j=0; j < n_across; j++) {
            // Point in physical space
            Vector2d pos = {r_pc[j], z_pc[i]};
            // Point in (z, Psi) space
            Vector2d pos_psi = {interp_rPsi.interpolated_value(pos), z_pc[i]};

            // Find angle of streamline to z-axis
            Vector2d gradPsi = interp_rPsi.gradient(pos, interp_Psi);
            gradPsi.normalize();
            double sinz = abs(gradPsi[1]);
            angle[i][j] = asin(sinz);

            // Find interpolated in (z, Psi) space values
            B_p[i][j] = interp_B_p.interpolated_value(pos_psi);
            B_phi[i][j] = interp_B_phi.interpolated_value(pos_psi);
            N[i][j] = interp_N.interpolated_value(pos_psi);
            Gamma[i][j] = interp_Gamma.interpolated_value(pos_psi);
            sigma[i][j] = interp_sigma.interpolated_value(pos_psi);
            theta[i][j] = interp_theta.interpolated_value(pos_psi);
            beta_phi[i][j] = interp_beta_phi.interpolated_value(pos_psi);
            jsq[i][j] = interp_jsq.interpolated_value(pos_psi);
            jsq_phi[i][j] = interp_jsq_phi.interpolated_value(pos_psi);
            Psi[i][j] = interp_Psi.interpolated_value(pos_psi);
        }
    }

    std::fstream fs;

    fs.open(mhd_run_name + "_angle_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, angle);
        fs.close();
    }


    fs.open(mhd_run_name + "_B_p_psi_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, B_p);
        fs.close();
    }

    fs.open(mhd_run_name + "_B_phi_psi_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, B_phi);
        fs.close();
    }

    fs.open(mhd_run_name + "_N_psi_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, N);
        fs.close();
    }

    fs.open(mhd_run_name + "_Gamma_psi_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, Gamma);
        fs.close();
    }

    fs.open(mhd_run_name + "_sigma_psi_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, sigma);
        fs.close();
    }

    fs.open(mhd_run_name + "_theta_psi_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, theta);
        fs.close();
    }

    fs.open(mhd_run_name + "_beta_phi_psi_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, beta_phi);
        fs.close();
    }

    fs.open(mhd_run_name + "_jsq_psi_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, jsq);
        fs.close();
    }

    fs.open(mhd_run_name + "_jsq_phi_psi_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, jsq_phi);
        fs.close();
    }

    fs.open(mhd_run_name + "_Psi_psi_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, Psi);
        fs.close();
    }
}


void check_psi_interpolations(const std::string& mhd_run_name) {

    Delaunay_triangulation tr_B_p, tr_B_phi, tr_N, tr_Gamma, tr_beta_phi, tr_jsq, tr_Psi, tr_rPsi;
    // Triangulations for (z, Psi, value)
    create_triangulation(mhd_run_name + "_B_p_field_psi.txt", &tr_B_p);
    create_triangulation(mhd_run_name + "_B_phi_field_psi.txt", &tr_B_phi);
    create_triangulation(mhd_run_name + "_Psi_field_psi.txt", &tr_Psi);
    create_triangulation(mhd_run_name + "_n_plasma_field_psi.txt", &tr_N);
    create_triangulation(mhd_run_name + "_Gamma_field_psi.txt", &tr_Gamma);
    create_triangulation(mhd_run_name + "_beta_phi_field_psi.txt", &tr_beta_phi);
    create_triangulation(mhd_run_name + "_jsq_plasma_field_psi.txt", &tr_jsq);
    // Triangulation for (z, r, Psi)
    create_triangulation(mhd_run_name + "_Psi_field.txt", &tr_rPsi);

    // Interpolates in (z, Psi) coordinates
    SimulationInterpolater interp_B_p(&tr_B_p);
    SimulationInterpolater interp_B_phi(&tr_B_phi);
    SimulationInterpolater interp_N(&tr_N);
    SimulationInterpolater interp_Gamma(&tr_Gamma, 1.0);
    SimulationInterpolater interp_beta_phi(&tr_beta_phi);
    SimulationInterpolater interp_jsq(&tr_jsq);
    SimulationInterpolater interp_Psi(&tr_Psi, 1.0);

    // Interpolate Psi(z, r): given (z, r) returns (Psi)
    SimulationInterpolater interp_rPsi(&tr_rPsi, 1.0);
    // Points where to find values
    size_t n_along = 5000;
    size_t n_across = 500;
    auto z_pc = linspace(0.0, 10.0, 5001);
    auto r_pc = linspace(0.0, 5.0, 501);
    auto psi = linspace(0.0, 1.0, 501);

    std::vector<std::vector<double>> B_p, B_phi, N, Gamma, beta_phi, jsq, Psi, angle;
    B_p.resize(n_along);
    B_phi.resize(n_along);
    N.resize(n_along);
    Gamma.resize(n_along);
    beta_phi.resize(n_along);
    jsq.resize(n_along);
    Psi.resize(n_along);
    angle.resize(n_along);
    for(size_t i=0; i < n_along; i++) {
        B_p[i].resize(n_across);
        B_phi[i].resize(n_across);
        N[i].resize(n_across);
        Gamma[i].resize(n_across);
        beta_phi[i].resize(n_across);
        jsq[i].resize(n_across);
        Psi[i].resize(n_across);
        angle[i].resize(n_across);
        for(size_t j=0; j < n_across; j++) {
            // Point in physical space
            Vector2d pos = {r_pc[j],z_pc[i]};
            // Point in (z, Psi) space
            Vector2d pos_psi = {interp_rPsi.interpolated_value(pos), z_pc[i]};

            // Find angle of streamline to z-axis
            Vector2d gradPsi = interp_rPsi.gradient(pos, interp_Psi);
            gradPsi.normalize();
            double sinz = abs(gradPsi[1]);
            angle[i][j] = asin(sinz);

            // Find interpolated in (z, Psi) space values
            B_p[i][j] = interp_B_p.interpolated_value(pos_psi);
            B_phi[i][j] = interp_B_phi.interpolated_value(pos_psi);
            N[i][j] = interp_N.interpolated_value(pos_psi);
            Gamma[i][j] = interp_Gamma.interpolated_value(pos_psi);
            beta_phi[i][j] = interp_beta_phi.interpolated_value(pos_psi);
            jsq[i][j] = interp_jsq.interpolated_value(pos_psi);
            Psi[i][j] = interp_Psi.interpolated_value(pos_psi);
        }
    }

    std::fstream fs;

    fs.open(mhd_run_name + "_angle_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, angle);
        fs.close();
    }


    fs.open(mhd_run_name + "_B_p_psi_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, B_p);
        fs.close();
    }

    fs.open(mhd_run_name + "_B_phi_psi_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, B_phi);
        fs.close();
    }

    fs.open(mhd_run_name + "_N_psi_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, N);
        fs.close();
    }

    fs.open(mhd_run_name + "_Gamma_psi_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, Gamma);
        fs.close();
    }

    fs.open(mhd_run_name + "_beta_phi_psi_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, beta_phi);
        fs.close();
    }

    fs.open(mhd_run_name + "_jsq_psi_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, jsq);
        fs.close();
    }

    fs.open(mhd_run_name + "_Psi_psi_interpolated.txt", std::ios::out | std::ios::app);
    if (fs.is_open()) {
        write_2dvector(fs, Psi);
        fs.close();
    }
}

std::vector<double> run_on_simulations(const std::string& mhd_run_name,
                                       double n_scale_nt, double n_scale_border, double n_scale_axis,
                                       double gamma_min, bool anisotropic_s,
                                       const std::string& particles_heating_model,
                                       double Psi_c_border, double sigma_Psi_border, double sigma_Psi_axis,
                                       double relerr) {

    auto t1 = Clock::now();
    std::clock_t start;
    start = std::clock();

    double los_angle = 17.0*M_PI/180.0;
    double redshift = 0.00436;

    // Observed frequencies in GHz
    //std::vector<double> nu_observed_ghz{15.4, 12.1, 8.1};
//    std::vector<double> nu_observed_ghz{15.4, 8.1};
//    std::vector<double> nu_observed_ghz{86.0};
    std::vector<double> nu_observed_ghz{15.4};
    std::vector<double> total_fluxes;
    // Frequencies in the BH frame in Hz
    std::vector<double> nu_bh;
    for(auto nu_obs_ghz : nu_observed_ghz) {
        nu_bh.push_back(nu_obs_ghz*1E+09*(1+redshift));
    }

    // Setting geometry using simulations ==============================================================================
    vector< vector<double> > all_points;
    read_from_txt(mhd_run_name + "_Psi_field.txt", all_points);
    size_t nrows = all_points.size();

    std::vector<Point_3> points;
    int n_circle = 36;
    std::cout << "Reading geometry file with #rows = " << nrows << std::endl;
    for (size_t i=0; i<nrows; i++) {
        // r & z in file are already in pc and we need in pc
        double z = all_points[i][0];
        double r_p = all_points[i][1];
//        std::cout << "z[pc] = " << z << ", r_p[pc] = " << r_p << "\n";
        for (int j=0; j<n_circle; j++) {
            double x = r_p*sin(j*2*M_PI/n_circle);
            double y = r_p*cos(j*2*M_PI/n_circle);
            double length_ = sqrt(x*x + y*y + z*z);
            points.emplace_back(Point_3(x, y, z));
        }
    }

    Polyhedron P;
    // FIXME: The problem is here
    CGAL::convex_hull_3(points.begin(), points.end(), P);
    std::cout << "Debug"<< "\n";

    Tree tree(faces(P).first, faces(P).second, P);
    SimulationGeometry geometry(&tree);


    // Setting interpolation of Psi
    Delaunay_triangulation tr_rPsi;
    // Triangulation for (z, r, Psi)
    create_triangulation(mhd_run_name + "_Psi_field.txt", &tr_rPsi);
    // Interpolate Psi(z, r): given (z, r) returns (Psi)
    SimulationInterpolater interp_rPsi(&tr_rPsi, 1.0);


    Delaunay_triangulation tr_Psi;
    // Triangulation for (Psi, r, Psi)
    create_triangulation(mhd_run_name + "_Psi_field_psi.txt", &tr_Psi);
    // Interpolate Psi(z, r): given (z, r) returns (Psi)
    SimulationInterpolater interp_Psi(&tr_Psi, 1.0);


    // Poloidal angle (in rad)
    Delaunay_triangulation tr_poloidal_angle;
    create_triangulation(mhd_run_name + "_theta_field_psi.txt", &tr_poloidal_angle);


    // Setting B-Field using simulations ===============================================================================
    Delaunay_triangulation tr_p;
    Delaunay_triangulation tr_fi;
    create_triangulation(mhd_run_name + "_B_p_field_psi.txt", &tr_p);
    create_triangulation(mhd_run_name + "_B_phi_field_psi.txt", &tr_fi);
    SimulationBField bfield(&tr_rPsi, &tr_Psi, &tr_poloidal_angle, &tr_p, &tr_fi, false);
//    std::vector<VectorBField*> vbfields;
//    vbfields.push_back(&bfield);


    // Setting N-field using simulations ===============================================================================
//    Delaunay_triangulation tr_nt;
    Delaunay_triangulation tr_ncold;
    Delaunay_triangulation tr_Bsq;
    Delaunay_triangulation tr_jsq;
    Delaunay_triangulation tr_sigma;

    create_triangulation(mhd_run_name + "_n_plasma_field_psi.txt", &tr_ncold);
    create_triangulation(mhd_run_name + "_jsq_plasma_field_psi.txt", &tr_jsq);
    create_triangulation(mhd_run_name + "_Bsq_plasma_field_psi.txt", &tr_Bsq);
    create_triangulation(mhd_run_name + "_sigma_field_psi.txt", &tr_sigma);


//    if(particles_heating_model == "n"){
//        // Number of emitting particles ~ number of cold particles
//        // In this case n_scale_nt must be in [0, 1]
//        create_triangulation(mhd_run_name + "_n_plasma_field_psi.txt", &tr_nt);
//    }
//    else if(particles_heating_model == "jsq") {
//        // Number of emitting particles ~ j_plasma^2
//        create_triangulation(mhd_run_name + "_jsq_plasma_field_psi.txt", &tr_nt);
//    }
//    else if(particles_heating_model == "bsq") {
//        // Number of emitting particles ~ B_plasma^2
//        // In this case n_scale_nt must be in [0, 1].
//        // u_e = n_nt * mc^2 * (s-1)/(s-2) * gamma_min = B^2/(8pi)
//        double s = 2.5;
//        // scale coefficient from B^2 to n_nt
//        double n_nt_Bsq = (s - 2)/(s - 1)/(8*M_PI*m_e*c*c*gamma_min);
//        n_scale_nt *= n_nt_Bsq;
//        create_triangulation(mhd_run_name + "_Bsq_plasma_field_psi.txt", &tr_nt);
//    } else {
//        std::vector<string> implemented_heating_model_types{"bsq", "n", "jsq"};
//        std::cout << "Implemented particles heating models are: " << "\n";
//        for(const auto& s: implemented_heating_model_types){
//            std::cout << s << "\n";
//        }
//        throw NotImplmentedParticlesHeating(particles_heating_model);
//    }

    // Arbitrary heating model w/o any constrains on the NT particles number density
//    SimulationNField nfield(&tr_nt, true, 2.5, gamma_min, anisotropic_s, n_scale_nt);

    // Heating model with constrain that number density of NT particles can not exceed some fraction of number density
    // of the cold particles
//    ConstrainedSimulationNField cnfield(&tr_nt, &tr_nt, true, 2.5, gamma_min, anisotropic_s, n_scale_nt);

    // Heating model with constrain that number density of NT particles can not exceed some fraction of number density
    // of the cold particles AND heating efficiency is modulated by local magnetization as in Broderick+2010
    double max_frac_cold = 0.1;
    double s = 2.5;
    std::string constrain_type;
//    constrain_type = "sigma";
    constrain_type = "none";
    std::cout << "Particles heating supression : " << constrain_type << "\n";
    // psi_mean - psi_width > 1 => no sheath!
//    ConstrainedBetaSimulationNField nfield(&tr_ncold, &tr_Bsq, &tr_jsq, &tr_sigma,
//                                           particles_heating_model, constrain_type, true,
//                                           s, gamma_min, anisotropic_s, n_scale_nt, max_frac_cold, n_scale_border,
//                                           0.5, 0.01);
    ByHandSimulationNField nfield(&tr_ncold, true, s, gamma_min, anisotropic_s,
                                  max_frac_cold, n_scale_border, Psi_c_border, sigma_Psi_border, sigma_Psi_axis,
                                  n_scale_axis, n_scale_nt);

    // Setting V-field using simulations ===============================================================================
    Delaunay_triangulation tr_Gamma;
    Delaunay_triangulation tr_beta_phi;
    create_triangulation(mhd_run_name + "_Gamma_field_psi.txt", &tr_Gamma);
    create_triangulation(mhd_run_name + "_beta_phi_field_psi.txt", &tr_beta_phi);
    SimulationVField vfield(&tr_rPsi, &tr_Psi, &tr_poloidal_angle, &tr_Gamma, &tr_beta_phi);

    Jet bkjet(&geometry, &interp_rPsi, &vfield, &bfield, &nfield);

    // FIXME: Put inside frequency loop for dep. on frequency
    // FIXME: Make coarser grid for scale estimation
    // Setting parameters of pixels and image ==========================================================================
    // These are OK for uniform pixel
//    int number_of_pixels_along = 600;
//    int number_of_pixels_across = 200;
//    double pixel_size_mas_start = 0.025;
//    double pixel_size_mas_stop = 0.025;
    // From 0.001 pc/pixel - that is for z=0.02 pc
    // Non-uniform pixel from ``pixel_size_mas_start`` (near BH) to ``pixel_size_mas_stop`` (image edges)
    int number_of_pixels_along = 800;
    int number_of_pixels_across = 150;
    double pixel_size_mas_start = 0.025;
    double pixel_size_mas_stop = 0.05;
    // 86 Ghz
//    int number_of_pixels_along = 600;
//    int number_of_pixels_across = 600;
//    double pixel_size_mas_start = 0.003;
//    double pixel_size_mas_stop = 0.003;

    auto image_size = std::make_pair(number_of_pixels_across, number_of_pixels_along);
    auto pc_in_mas = mas_to_pc(redshift);
    std::cout << "pc_in_mas " << pc_in_mas << std::endl;
    // Log10 of pixel size in cm
    auto lg_pixel_size_start = log10(pixel_size_mas_start*pc_in_mas*pc);
    auto lg_pixel_size_stop = log10(pixel_size_mas_stop*pc_in_mas*pc);

    std::cout << "Setting pixel size (pc) from " << pow(10.0, lg_pixel_size_start)/pc << " to " << pow(10.0, lg_pixel_size_stop)/pc << std::endl;
    for(auto jet_side : {true, false}) {

        // Ignore CJ
        if(jet_side == false) {
            continue;
        }

        ImagePlane imagePlane(image_size, lg_pixel_size_start, lg_pixel_size_stop, los_angle, jet_side);
        // Array of pixel sizes in cm
        auto pixel_sizes = imagePlane.getPixelSizes();
        // Array of pixel solid angles in rad*rad
        std::vector<std::vector<double>> pixel_solid_angles;
        pixel_solid_angles.resize(pixel_sizes.size());

        for(unsigned long i=0; i < pixel_sizes.size(); i++) {
            pixel_solid_angles[i].resize(pixel_sizes[0].size());
            for(unsigned long j=0; j < pixel_sizes[0].size(); j++) {
                // Divide by ``pc_in_mas*pc`` to bring ``cm`` to ``mas`` at source redshift
                pixel_solid_angles[i][j] = (pixel_sizes[i][j]/(pc_in_mas*pc))*(pixel_sizes[i][j]/(pc_in_mas*pc))*mas_to_rad*mas_to_rad;
            }
        }

        // Array of scale factors. Divide resulting image on this to obtain flux density in Jy. Accounts for cosmological
        // scaling of intensity
        std::vector<std::vector<double>> scales;
        scales.resize(pixel_sizes.size());
        for(unsigned long i=0; i < pixel_sizes.size(); i++) {
            scales[i].resize(pixel_sizes[0].size());
            for(unsigned long j=0; j < pixel_sizes[0].size(); j++) {
                scales[i][j] = 1E-23*(1.+redshift)*(1.+redshift)*(1.+redshift)/pixel_solid_angles[i][j];
//                scales[i][j] = (1.+redshift)*(1.+redshift)*(1.+redshift)/pixel_solid_angles[i][j];
            }
        }

        Observation observation(&bkjet, &imagePlane);

        // FIXME: Put out of frequency loop - these do not depend on frequency
        // Transfer-specific parameters ================================================================================
        double tau_max = 10;
        double dt_max_pc = 0.01;
        double dt_max = pc*dt_max_pc;
        double tau_min_log10 = -20.0;
        double tau_min = pow(10.,tau_min_log10);
        int n_ = 100;

        string polarization = "I";
//        string polarization = "speed";
//        string polarization = "full";

        for(int i_nu=0; i_nu < nu_observed_ghz.size(); i_nu++) {

            if(jet_side) {
                std::cout << "Running transfer for frequency " << nu_observed_ghz[i_nu] << " GHz for approaching jet" << std::endl;
            } else {
                std::cout << "Running transfer for frequency " << nu_observed_ghz[i_nu] << " GHz for counter-jet" << std::endl;
            }
            observation.run(n_, tau_max, dt_max, tau_min, nu_bh[i_nu], polarization, relerr);
            string value = "tau";
            auto image_tau = observation.getImage(value);

            value = "I";
            auto image_i = observation.getImage(value);
            double total_flux = 0.0;
            for (unsigned long int i = 0; i < image_i.size(); ++i) {
                for (unsigned long int j = 0; j < image_i[i].size(); ++j) {
                    image_i[i][j] = image_i[i][j]/scales[i][j];
                    total_flux += image_i[i][j];
                }
            }

            if(jet_side == true){
                total_fluxes.push_back(total_flux);
            }

            value = "l";
            auto image_l = observation.getImage(value);

            std::fstream fs;
            // Remove trailing zeros: https://stackoverflow.com/a/46424921
            std::ostringstream oss;
            oss << std::setprecision(8) << std::noshowpoint << nu_observed_ghz[i_nu];
            std::string freq_name = oss.str();

            std::string prefix = "_psi_" + std::to_string(Psi_c_border) + "_dpsi_" + std::to_string(sigma_Psi_border);

            std::string file_tau, file_tau_fr, file_i, file_beta, file_q, file_u, file_v, file_l;
            if(jet_side) {
                file_tau = mhd_run_name + prefix + "_jet_image_tau_" + freq_name + ".txt";
                file_tau_fr = mhd_run_name + prefix + "_jet_image_taufr_" + freq_name + ".txt";
                file_i = mhd_run_name + prefix + "_jet_image_i_" + freq_name + ".txt";
                file_beta = mhd_run_name + prefix + "_jet_image_beta_app_" + freq_name + ".txt";
                file_q = mhd_run_name + prefix + "_jet_image_q_" + freq_name + ".txt";
                file_u = mhd_run_name + prefix + "_jet_image_u_" + freq_name + ".txt";
                file_v = mhd_run_name + prefix + "_jet_image_v_" + freq_name + ".txt";
                file_l = mhd_run_name + prefix + "_jet_image_l_" + freq_name + ".txt";
            } else {
                file_tau = mhd_run_name + "_cjet_image_tau_" + freq_name + ".txt";
                file_tau_fr = mhd_run_name + "_cjet_image_taufr_" + freq_name + ".txt";
                file_i = mhd_run_name + "_cjet_image_i_" + freq_name + ".txt";
                file_beta = mhd_run_name + "_cjet_image_beta_app_" + freq_name + ".txt";
                file_q = mhd_run_name + "_cjet_image_q_" + freq_name + ".txt";
                file_u = mhd_run_name + "_cjet_image_u_" + freq_name + ".txt";
                file_v = mhd_run_name + "_cjet_image_v_" + freq_name + ".txt";
                file_l = mhd_run_name + "_cjet_image_l_" + freq_name + ".txt";
            }

            // Remove old file
            std::remove(file_i.c_str());
            std::remove(file_beta.c_str());
            std::remove(file_q.c_str());
            std::remove(file_u.c_str());
            std::remove(file_v.c_str());
            std::remove(file_l.c_str());
            std::remove(file_tau.c_str());
            std::remove(file_tau_fr.c_str());

            fs.open(file_tau, std::ios::out | std::ios::app);
            if (fs.is_open()) {
                write_2dvector(fs, image_tau);
                fs.close();
            }

            fs.open(file_i, std::ios::out | std::ios::app);
            if (fs.is_open()) {
                write_2dvector(fs, image_i);
                fs.close();
            }

            fs.open(file_l, std::ios::out | std::ios::app);
            if (fs.is_open()) {
                write_2dvector(fs, image_l, pc);
                fs.close();
            }

            if (polarization == "speed") {
                value = "SPEED";
                auto image_beta_app = observation.getImage(value);
                // Scale emission weighted beta_app on total intensity in given pixel
                for (unsigned long int i = 0; i < image_beta_app.size(); ++i) {
                    for (unsigned long int j = 0; j < image_beta_app[i].size(); ++j) {
                        // ``image_i`` was already devided by scales. Need to recover the original intensities.
                        image_beta_app[i][j] = image_beta_app[i][j]/(image_i[i][j]*scales[i][j]);
                    }
                }
                fs.open(file_beta, std::ios::out | std::ios::app);
                if (fs.is_open()) {
                    write_2dvector(fs, image_beta_app);
                    fs.close();
                }
            }

            if (polarization == "full") {
                value = "Q";
                auto image_q = observation.getImage(value);
                for (unsigned long int i = 0; i < image_q.size(); ++i) {
                    for (unsigned long int j = 0; j < image_q[i].size(); ++j) {
                        image_q[i][j] = image_q[i][j]/scales[i][j];
                    }
                }

                value = "U";
                auto image_u = observation.getImage(value);
                for (unsigned long int i = 0; i < image_u.size(); ++i) {
                    for (unsigned long int j = 0; j < image_u[i].size(); ++j) {
                        image_u[i][j] = image_u[i][j]/scales[i][j];
                    }
                }

                value = "V";
                auto image_v = observation.getImage(value);
                for (unsigned long int i = 0; i < image_v.size(); ++i) {
                    for (unsigned long int j = 0; j < image_v[i].size(); ++j) {
                        image_v[i][j] = image_v[i][j]/scales[i][j];
                    }
                }

                value = "tau_fr";
                auto image_tau_fr = observation.getImage(value);

                fs.open(file_tau_fr, std::ios::out | std::ios::app);
                if (fs.is_open()) {
                    write_2dvector(fs, image_tau_fr);
                    fs.close();
                }

                fs.open(file_q, std::ios::out | std::ios::app);
                if (fs.is_open()) {
                    write_2dvector(fs, image_q);
                    fs.close();
                }

                fs.open(file_u, std::ios::out | std::ios::app);
                if (fs.is_open()) {
                    write_2dvector(fs, image_u);
                    fs.close();
                }

                fs.open(file_v, std::ios::out | std::ios::app);
                if (fs.is_open()) {
                    write_2dvector(fs, image_v);
                    fs.close();
                }
            }
        }
    }

    std::cout << "CPU Time: "
              << (std::clock() - start) / (double) (CLOCKS_PER_SEC)
              << " s" << std::endl;
    auto t2 = Clock::now();
    std::cout << "User time: "
              << std::chrono::duration_cast<std::chrono::seconds>(
                  t2 - t1).count()
              << " s" << std::endl;

    return total_fluxes;
}


// To run in parallel when fil params_3ridges.txt has 3 parameter sets:
// parallel --files --results result_{1}_AmpAxis_{3}_AmpEdges_{4}_PsiEdges_{5}_dPsiEdges_{6}_dPsiAxis_{7}_relerr_{8} --joblog log --jobs 3 -a params_3ridges.txt -n 1 -m --colsep ' ' "./mhd_transfer"
int main(int argc, char *argv[]) {

    bool anisotropic_s = false;
    std::string particles_heating_model = "byhand";
    double gamma_min = 100.0;

    std::vector<string> implemented_heating_model_types{"bsq", "n", "jsq", "byhand"};
    std::vector<double> total_fluxes;

//    if(argc != 9){
//        std::cout << argc << "\n";
//        std::cout << "Supply MHD code, NT uniform density scale factor (0 < f [< 1]), NT axis scale factor (0 < f [<1]),"
//                     " NT border density scale factor (0 < f [< 1]),"
//                     " center Psi for border, width Psi for border, width Psi for axis, rel.error\n" << "\n";
//        return 1;
//    }
//    else {
//        std::cout << "Doing radiation transport for MHD code " << argv[1] << "\n";
//
//        double n_scale_nt = atof(argv[2]);
//        std::cout << "scaling factor for NT particles uniform density = " << argv[2] << "\n";
//
//        double n_scale_axis = atof(argv[3]);
//        std::cout << "scaling factor for NT particles density at the axis = " << argv[3] << "\n";
//
//        double n_scale_border = atof(argv[4]);
//        std::cout << "scaling factor for NT particles density at the border = " << argv[4] << "\n";
//
//        double Psi_c = atof(argv[5]);
//        std::cout << "Psi_c for border = " << argv[5] << "\n";
//
//        double sigma_Psi = atof(argv[6]);
//        std::cout << "sigma_Psi for border = " << argv[6] << "\n";
//
//        double sigma_Psi_axis = atof(argv[7]);
//        std::cout << "sigma_Psi for axis = " << argv[7] << "\n";
//
//        double relerr = atof(argv[8]);
//        std::cout << "relerr = " << argv[8] << "\n";
//
////        double gamma_min = atof(argv[4]);
////        std::cout << "gamma_min = " << argv[4] << "\n";
//
////        bool anisotropic_s;
////        std::istringstream is(argv[5]);
////        is >> std::boolalpha >> anisotropic_s;
////        std::cout << "Anisotropic s = " << argv[5] << "\n";
//
////        std::string particles_heating_model = argv[5];
////        std::cout << "Particles heating model = " << argv[5] << "\n";
//        // Check that particles heating model is implemented. This check is also done inside ``run_on_simulations``
//        if (std::find(implemented_heating_model_types.begin(),
//                      implemented_heating_model_types.end(),
//                      particles_heating_model) == implemented_heating_model_types.end())
//        {
//            throw NotImplmentedParticlesHeating(particles_heating_model);
//        }
//        total_fluxes = run_on_simulations(argv[1], n_scale_nt, n_scale_border, n_scale_axis, gamma_min,
//                                          anisotropic_s, particles_heating_model,
//                                          Psi_c, sigma_Psi, sigma_Psi_axis,
//                                          relerr);
//    }
//    for(auto total_flux: total_fluxes){
//        std::cout << "Total flux [Jy] = " << total_flux << "\n";
//    }
//
//
    // From IDE
    total_fluxes = run_on_simulations("m1s10g2b123.971372r0.000369", 0.0, 4.0, 2.0,
                                      gamma_min, anisotropic_s, particles_heating_model,
                                      1.0, 0.015, 0.03, 1e-05);
    for(auto total_flux: total_fluxes){
        std::cout << "Total flux [Jy] = " << total_flux << "\n";
    }

//    check_psi_interpolations_Lena("m1s10g2b123.971372r0.000369");
//    check_psi_interpolations_Lena("m2s10g2b44.614955r0.000595");
//    check_psi_interpolations("eta");
    return 0;
}