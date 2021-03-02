#ifndef MHD_TRANSFER_BFIELD_H
#define MHD_TRANSFER_BFIELD_H

#include <Eigen/Eigen>
#include "SimulationInterpolater.h"
#include "Geometry.h"

using Eigen::Vector3d;


// B-field with vector values, e.g. ordered component or ordered component with cells, with specified fraction of the
// completely tangled component. Vector component can be specified in any frame (plasma or lab). Tangled component is
// specified only in plasma frame as some fraction of the vector component.
class VectorBField {
    public:
        virtual Vector3d _bf(const Vector3d &point, double psi) const = 0 ;
        Vector3d bf(const Vector3d &point, double psi) const ;
        // B-field in plasma (comoving) frame. Needed for calculation of transfer coefficients
        Vector3d bf_plasma_frame(const Vector3d &point, double psi, Vector3d &v) const;
        // Tangled B-field component in plasma (comoving) frame. Needed for calculation of transfer coefficients
        double bf_tangled_plasma_frame(const Vector3d &point, double psi, Vector3d &v) const;
        // Unit vector of B-field in laboratory (observer) frame. Needed for calculation of polarization swing.
        Vector3d bhat_lab_frame(const Vector3d &point, double psi, Vector3d &v) const;
        double get_tangled_fraction() const {
            return tangled_fraction_;
        };

    protected:
        VectorBField(bool in_plasma_frame, double tangled_fraction, Geometry* geometry_out=nullptr, Geometry* geometry_in=nullptr);
        bool in_plasma_frame_;
        double tangled_fraction_;
        Geometry* geometry_in_;
        Geometry* geometry_out_;

};


// Special class for Elena's calculations of axisymmetric flows.
class SimulationBField : public VectorBField {
    public:
        SimulationBField(Delaunay_triangulation *tr_psi, Delaunay_triangulation *tr_p, Delaunay_triangulation *tr_fi,
                         bool in_plasma_frame, double tangled_fraction=0.0);
        Vector3d _bf(const Vector3d &point, double psi) const override ;

    private:
        SimulationInterpolater interp_p_;
        SimulationInterpolater interp_fi_;
        // For calculating direction of the poloidal component we need Psi(r_p, z)
        SimulationInterpolater interp_psi_;
};


#endif //MHD_TRANSFER_BFIELD_H
