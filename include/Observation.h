#ifndef MHD_TRANSFER_OBSERVATION_H
#define MHD_TRANSFER_OBSERVATION_H

#include <string>
#include "Jet.h"
#include "ImagePlane.h"


class Observation {
    public:
        Observation(Jet* newjet, ImagePlane* imagePlane);
        void run(int n, double tau_max, double dt_max, double tau_min, double nu, string polarization);
        void observe_single_pixel(Ray& ray, Pixel& pixel, double tau_min, double tau_max, int n, double dt_max,
                                  double nu, string polarization);
        std::vector<std::vector<double>> getImage(string value);
        std::pair<unsigned long int,unsigned long int> getImageSize();
    private:
        Jet* jet;
        ImagePlane* imagePlane;

		std::pair<double, double> integrate_tau_adaptive(std::list<Intersection>& list_intersect, Vector3d ray_direction,
		                                                 const double nu, double tau_max, int n, double dt_max);

		void integrate_i_adaptive(std::list<Intersection>& list_intersect, const Vector3d& ray_direction, const double nu,
		                          int n, double& background_I, double dt_max);

        void integrate_full_stokes(std::list<Intersection>& list_intersect, const Vector3d& ray_direction, const double nu,
                                   int n, std::vector<double>& background, double dt_max);

		void integrate_full_stokes_adaptive(std::list<Intersection>& list_intersect, const Vector3d& ray_direction,
		                                    const double nu, int n, std::vector<double>& background, double dt_max);

        void integrate_faraday_rotation_depth_adaptive(std::list<Intersection>& list_intersect, const Vector3d& ray_direction,
                                                       const double nu, int n, double& background, double dt_max);

};

#endif //MHD_TRANSFER_OBSERVATION_H
