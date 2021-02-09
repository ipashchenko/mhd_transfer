#ifndef MHD_TRANSFER_VFIELD_H
#define MHD_TRANSFER_VFIELD_H

#include <Eigen/Eigen>
#include "SimulationInterpolater.h"


using Eigen::Vector3d;

// TODO: Add ``gamma`` getter - Jet instances need it if particles density in plasma frame must be calculated
class VField {
    public:
        virtual Vector3d vf(const Vector3d &point) const = 0;
};


class SimulationVField: public VField {
    public:
        explicit SimulationVField(Delaunay_triangulation *tr_Gamma, Delaunay_triangulation *tr_beta_phi);
        Vector3d vf(const Vector3d& point) const override;

    private:
        SimulationInterpolater interp_Gamma_;
        SimulationInterpolater interp_beta_phi_;
};

#endif //MHD_TRANSFER_VFIELD_H