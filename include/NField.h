#ifndef MHD_TRANSFER_NFIELD_H
#define MHD_TRANSFER_NFIELD_H


#include <Eigen/Eigen>
#include "SimulationInterpolater.h"

using Eigen::Vector3d;


class NField {
    public:
        // This returns density
        virtual double _nf(const Vector3d &point) const = 0;
        // This wraps possible Geometry
        double nf(const Vector3d &point) const;
        double nf_plasma_frame(const Vector3d &point, double &gamma) const;

        virtual double k_i(double b, Vector3d &n_los, double nu, double n) const = 0;
        // Absorption coefficient for given vector of magnetic field ``b``, unit LOS
        // vector ``n_los`` and others measured in emission frame
        virtual double k_i(Vector3d &b, Vector3d &n_los, double nu, double n) const = 0;
        virtual double eta_i(double b, Vector3d &n_los, double nu, double n) const = 0;
        // Emission coefficient for given vector of magnetic field ``b``, unit LOS
        // vector ``n_los`` and others measured in emission frame
        virtual double eta_i(Vector3d &b, Vector3d &n_los, double nu, double n) const = 0;
        virtual double k_q(Vector3d &b, Vector3d &n_los, double nu, double n) const = 0;
        //virtual double k_q(double b, Vector3d &n_los, double nu, double n) const = 0;
        virtual double k_u(Vector3d &b, Vector3d &n_los, double nu, double n) const = 0;
        //virtual double k_u(double b, Vector3d &n_los, double nu, double n) const = 0;
        virtual double k_v(Vector3d &b, Vector3d &n_los, double nu, double n) const = 0;
        //virtual double k_v(double b, Vector3d &n_los, double nu, double n) const = 0;

        double k_F_c(Vector3d &b, Vector3d &n_los, double nu, double n) const;
        //double k_F_c(double b, Vector3d &n_los, double nu, double n) const;
        double k_C_c(Vector3d &b, Vector3d &n_los, double nu, double n) const;
        //double k_C_c(double b, Vector3d &n_los, double nu, double n) const;

        virtual double k_F(Vector3d &b, Vector3d &n_los, double nu, double n) const = 0;
        //virtual double k_F(double b, Vector3d &n_los, double nu, double n) const = 0;
        virtual double k_C(Vector3d &b, Vector3d &n_los, double nu, double n) const = 0;
        //virtual double k_C(double b, Vector3d &n_los, double nu, double n) const = 0;
        virtual double h_Q(Vector3d &b, Vector3d &n_los, double nu, double n) const = 0;
        //virtual double h_Q(double b, Vector3d &n_los, double nu, double n) const = 0;
        virtual double eta_q(Vector3d &b, Vector3d &n_los, double nu, double n) const = 0;
        //virtual double eta_q(double b, Vector3d &n_los, double nu, double n) const = 0;
        virtual double eta_u(Vector3d &b, Vector3d &n_los, double nu, double n) const = 0;
        //virtual double eta_u(double b, Vector3d &n_los, double nu, double n) const = 0;
        virtual double eta_v(Vector3d &b, Vector3d &n_los, double nu, double n) const = 0;
        //virtual double eta_v(double b, Vector3d &n_los, double nu, double n) const = 0;

    protected:
        explicit NField(bool in_plasma_frame, Geometry* geometry = nullptr);
        bool in_plasma_frame_;
        Geometry* geometry_;
};


class PowerLawNField : public NField {
    public:
        PowerLawNField(bool in_plasma_frame, double s, double gamma_min, Geometry* geometry = nullptr,
                       std::string plasma="normal", bool changing_s="false", double ds=0.0);
        double get_s(const Vector3d &point) const;
        double get_s(Vector3d &b, Vector3d &n_los) const;

        double k_i(double b, Vector3d &n_los, double nu, double n) const override;
        double k_i(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        double eta_i(double b, Vector3d &n_los, double nu, double n) const override;
        double eta_i(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        double k_q(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        //double k_q(double b, Vector3d &n_los, double nu, double n) const override;
        double k_u(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        //double k_u(double b, Vector3d &n_los, double nu, double n) const override;
        double k_v(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        //double k_v(double b, Vector3d &n_los, double nu, double n) const override;
        double k_F(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        //double k_F(double b, Vector3d &n_los, double nu, double n) const override;
        double k_C(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        //double k_C(double b, Vector3d &n_los, double nu, double n) const override;
        double h_Q(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        //double h_Q(double b, Vector3d &n_los, double nu, double n) const override;
        double eta_q(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        //double eta_q(double b, Vector3d &n_los, double nu, double n) const override;
        double eta_u(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        //double eta_u(double b, Vector3d &n_los, double nu, double n) const override;
        double eta_v(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        //double eta_v(double b, Vector3d &n_los, double nu, double n) const override;
    private:
        std::string plasma_;
        double s_;
        double gamma_min_;
        bool changing_s_;
        double ds_;
        double factor_ki_;
        double factor_ki_rnd_;
        double factor_kv_;
        double factor_kf_;
        double factor_etai_;
        double factor_etai_rnd_;
        double factor_etav_;
};


class ThermalNField : public NField {
    public:
        ThermalNField(Geometry* geometry);
        double k_i(double b, Vector3d &n_los, double nu, double n) const override;
        double k_i(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        double eta_i(double b, Vector3d &n_los, double nu, double n) const override;
        double eta_i(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        double k_q(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        //double k_q(double b, Vector3d &n_los, double nu, double n) const override;
        double k_u(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        //double k_u(double b, Vector3d &n_los, double nu, double n) const override;
        double k_v(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        //double k_v(double b, Vector3d &n_los, double nu, double n) const override;
        double k_F(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        //double k_F(double b, Vector3d &n_los, double nu, double n) const override;
        double k_C(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        //double k_C(double b, Vector3d &n_los, double nu, double n) const override;
        double h_Q(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        //double h_Q(double b, Vector3d &n_los, double nu, double n) const override;
        double eta_q(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        //double eta_q(double b, Vector3d &n_los, double nu, double n) const override;
        double eta_u(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        //double eta_u(double b, Vector3d &n_los, double nu, double n) const override;
        double eta_v(Vector3d &b, Vector3d &n_los, double nu, double n) const override;
        //double eta_v(double b, Vector3d &n_los, double nu, double n) const override;
};





class SimulationNField: public PowerLawNField {
    public:
        SimulationNField(Delaunay_triangulation *tr, bool in_plasma_frame, double s, double gamma_min,
                         bool changing_s=false, double scale_factor=1.0);
        double _nf(const Vector3d &point) const override;

    private:
        double scale_factor_;
        SimulationInterpolater interp_;
};


#endif //MHD_TRANSFER_NFIELD_H
