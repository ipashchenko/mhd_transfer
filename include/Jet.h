#ifndef MHD_TRANSFER_JET_H
#define MHD_TRANSFER_JET_H


#include "Geometry.h"
#include "VField.h"
#include "BField.h"
#include "NField.h"
#include "utils.h"


class Jet {
    public:
        // Need VectorBField for simulations output or polarization
        Jet(BaseGeometry* geo, VField* vfield, VectorBField* bField, NField* nField);
        // Need ScalarBField for analytical Blandford-Konigle
        //Jet(BaseGeometry* geo, VField* vfield, ScalarBField* bField, NField* nField);

        // Get all coefficients
        std::tuple<double, double, double, double, double, double, double, double, double, double, double> get_transport_coefficients(Vector3d &point, Vector3d &n_los, double nu);

        // Absorption coefficient in ``point`` of the jet in the observer (lab)
        // frame. ``n`` is unit LOS vector in the observer frame.
        double getKI(Vector3d &point, Vector3d &n_los, double nu);

        // Emission coefficient in ``point`` of the jet in the observer (lab)
        // frame. ``n`` is unit LOS vector in the observer frame.
        double getEtaI(Vector3d &point, Vector3d &n_los, double nu);

        double getKF(Vector3d &point, Vector3d &n_los, double nu);

        std::list<Intersection> hit(Ray& ray);

        // Vector of the bulk motion speed (in cm/s) in the lab frame at point ``point``.
        Vector3d getV(const Vector3d& point);

        // Vector of the bulk motion speed (in c units) in the lab frame at point ``point``.
        const Vector3d getBeta(const Vector3d& point);

        // B-field in the lab frame at point ``point``.
        const Vector3d getB(const Vector3d& point);

        // Unit vector of B-field in the lab frame at point ``point``.
        const Vector3d getBhat(const Vector3d& point);


    private:
        BaseGeometry* geometry_;
        VField* vfield_;
        VectorBField* bfield_;
        NField* nfield_;
};



#endif //MHD_TRANSFER_JET_H
