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

using Eigen::Vector3d;
using Eigen::Matrix3Xd;
using std::vector;
using std::pair;
using namespace boost::numeric::odeint;
namespace ph = std::placeholders;

// For running w/o simulations
typedef std::chrono::high_resolution_clock Clock;


void run_on_simulations() {
    auto t1 = Clock::now();
    std::clock_t start;
    start = std::clock();

    double los_angle = 17.0*M_PI/180.0;
    double redshift = 0.00436;

    // Observed frequencies in GHz
    //std::vector<double> nu_observed_ghz{15.4, 12.1, 8.1};
    std::vector<double> nu_observed_ghz{15.4};
    // Frequencies in the BH frame in Hz
    std::vector<double> nu_bh;
    for(auto nu_obs_ghz : nu_observed_ghz) {
        nu_bh.push_back(nu_obs_ghz*1E+09*(1+redshift));
    }

    // Setting geometry using simulations ==============================================================================
    vector< vector<double> > all_points;
    read_from_txt("Gamma_field.txt", all_points);
    size_t nrows = all_points.size();

    std::vector<Point_3> points;
    int n_circle = 36;
    std::cout << "Reading geometry file with #rows = " << nrows << std::endl;
    for (size_t i=0; i<nrows; i++) {
        double z = all_points[i][0]/pc;
        double r_p = all_points[i][1]/pc;
        for (int j=0; j<n_circle; j++) {
            double x = r_p*sin(j*2*pi/n_circle);
            double y = r_p*cos(j*2*pi/n_circle);
            double length_ = sqrt(x*x + y*y + z*z);
            points.emplace_back(Point_3(x, y, z));
        }
    }

    Polyhedron P;
    CGAL::convex_hull_3(points.begin(), points.end(), P);
    Tree tree(faces(P).first, faces(P).second, P);
    SimulationGeometry geometry(&tree);

    // Setting B-Field using simulations ===============================================================================
    Delaunay_triangulation tr_p;
    Delaunay_triangulation tr_fi;
    create_triangulation("B_p_field.txt", &tr_p);
    create_triangulation("B_phi_field.txt", &tr_fi);
    SimulationBField bfield(&tr_p, &tr_fi, false);
//    std::vector<VectorBField*> vbfields;
//    vbfields.push_back(&bfield);



    // Setting N-field using simulations ===============================================================================
    Delaunay_triangulation tr_n;
    create_triangulation("n_plasma_field.txt", &tr_n);
    SimulationNField nfield(&tr_n, false, 2.5, 100.0);
//    std::vector<NField*> nfields;
//    nfields.push_back(&nfield);

    // Setting V-field using simulations ===============================================================================
    Delaunay_triangulation tr_v;
    create_triangulation("Gamma_field.txt", &tr_v);
    SimulationVField vfield(&tr_v);

    Jet bkjet(&geometry, &vfield, &bfield, &nfield);

    // FIXME: Put inside frequency loop for dep. on frequency
    // Setting parameters of pixels and image ==========================================================================
    int number_of_pixels_along = 300;
    int number_of_pixels_across = 80;
    // Non-uniform pixel from ``pixel_size_mas_start`` (near BH) to ``pixel_size_mas_stop`` (image edges)
    //double pixel_size_mas_start = 0.02;
    //double pixel_size_mas_stop = 0.5;
    double pixel_size_mas_start = 0.05;
    double pixel_size_mas_stop = 0.5;
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
            }
        }

        Observation observation(&bkjet, &imagePlane);

        // FIXME: Put out of frequency loop - these do not depend on frequency
        // Transfer-specific parameters ================================================================================
        double tau_max = 10;
        double dt_max_pc = 0.001;
        double dt_max = pc*dt_max_pc;
        double tau_min_log10 = -10.0;
        double tau_min = pow(10.,tau_min_log10);
        int n_ = 100;

        string polarization = "full";

        for(int i_nu=0; i_nu < nu_observed_ghz.size(); i_nu++) {

            if(jet_side) {
                std::cout << "Running transfer for frequency " << nu_observed_ghz[i_nu] << " GHz for approaching jet" << std::endl;
            } else {
                std::cout << "Running transfer for frequency " << nu_observed_ghz[i_nu] << " GHz for counter-jet" << std::endl;
            }

            observation.run(n_, tau_max, dt_max, tau_min, nu_bh[i_nu], polarization);
            string value = "tau";
            auto image_tau = observation.getImage(value);

            value = "I";
            auto image_i = observation.getImage(value);
            for (unsigned long int i = 0; i < image_i.size(); ++i) {
                for (unsigned long int j = 0; j < image_i[i].size(); ++j) {
                    image_i[i][j] = image_i[i][j]/scales[i][j];
                }
            }

            value = "l";
            auto image_l = observation.getImage(value);



            std::fstream fs;
            // Remove trailing zeros: https://stackoverflow.com/a/46424921
            std::ostringstream oss;
            oss << std::setprecision(8) << std::noshowpoint << nu_observed_ghz[i_nu];
            std::string freq_name = oss.str();

            std::string file_tau, file_tau_fr, file_i, file_q, file_u, file_v, file_l;
            if(jet_side) {
                file_tau = "jet_image_tau_" + freq_name + ".txt";
                file_tau_fr = "jet_image_taufr_" + freq_name + ".txt";
                file_i = "jet_image_i_" + freq_name + ".txt";
                file_q = "jet_image_q_" + freq_name + ".txt";
                file_u = "jet_image_u_" + freq_name + ".txt";
                file_v = "jet_image_v_" + freq_name + ".txt";
                file_l = "jet_image_l_" + freq_name + ".txt";
            } else {
                file_tau = "cjet_image_tau_" + freq_name + ".txt";
                file_tau_fr = "cjet_image_taufr_" + freq_name + ".txt";
                file_i = "cjet_image_i_" + freq_name + ".txt";
                file_q = "cjet_image_q_" + freq_name + ".txt";
                file_u = "cjet_image_u_" + freq_name + ".txt";
                file_v = "cjet_image_v_" + freq_name + ".txt";
                file_l = "cjet_image_l_" + freq_name + ".txt";
            }

            // Remove old file
            std::remove(file_i.c_str());
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
}


int main() {
    run_on_simulations();
    return 0;
}