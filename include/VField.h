#ifndef MHD_TRANSFER_VFIELD_H
#define MHD_TRANSFER_VFIELD_H

#include <Eigen/Eigen>
#include "SimulationInterpolater.h"


using Eigen::Vector3d;

// TODO: Add ``gamma`` getter - Jet instances need it if particles density in plasma frame must be calculated
class VField {
    public:
        // Psi needed for interpolation, but could be also helpful for analytical field
        virtual Vector3d vf(const Vector3d &point, double psi) const = 0;
};


class SimulationVField: public VField {
    public:
        explicit SimulationVField(Delaunay_triangulation *tr_rpsi,
                                  Delaunay_triangulation *tr_psi,
                                  Delaunay_triangulation *tr_poloidal_angle,
                                  Delaunay_triangulation *tr_Gamma,
                                  Delaunay_triangulation *tr_beta_phi);
        Vector3d vf(const Vector3d& point, double psi) const override;

    private:
        // Interpolators from (Psi, z)
        SimulationInterpolater interp_Gamma_;
        SimulationInterpolater interp_beta_phi_;
        // Already calculated poloidal angle interpolation
        SimulationInterpolater interp_polangle_;
        // For calculating direction of the poloidal component we need both Psi(r_p, z) and Psi(Psi, z)
        SimulationInterpolater interp_rpsi_;
        SimulationInterpolater interp_psi_;
};

#endif //MHD_TRANSFER_VFIELD_H
