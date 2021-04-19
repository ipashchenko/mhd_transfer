#ifndef MHD_TRANSFER_NFIELD_H
#define MHD_TRANSFER_NFIELD_H


#include <Eigen/Eigen>
#include "SimulationInterpolater.h"

using Eigen::Vector3d;


class NField {
    public:
        // This returns density
        virtual double _nf(const Vector3d &point, double psi) const = 0;
        // This wraps possible Geometry
        double nf(const Vector3d &point, double psi) const;
        double nf_plasma_frame(const Vector3d &point, double psi, double &gamma) const;

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

        double k_i(double b, Vector3d &n_los, double nu, double n_nt) const override;
        double k_i(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const override;
        double eta_i(double b, Vector3d &n_los, double nu, double n_nt) const override;
        double eta_i(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const override;
        double k_q(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const override;
        //double k_q(double b, Vector3d &n_los, double nu, double n_nt) const override;
        double k_u(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const override;
        //double k_u(double b, Vector3d &n_los, double nu, double n_nt) const override;
        double k_v(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const override;
        //double k_v(double b, Vector3d &n_los, double nu, double n_nt) const override;
        double k_F(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const override;
        //double k_F(double b, Vector3d &n_los, double nu, double n_nt) const override;
        double k_C(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const override;
        //double k_C(double b, Vector3d &n_los, double nu, double n_nt) const override;
        double h_Q(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const override;
        //double h_Q(double b, Vector3d &n_los, double nu, double n_nt) const override;
        double eta_q(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const override;
        //double eta_q(double b, Vector3d &n_los, double nu, double n_nt) const override;
        double eta_u(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const override;
        //double eta_u(double b, Vector3d &n_los, double nu, double n_nt) const override;
        double eta_v(Vector3d &b, Vector3d &n_los, double nu, double n_nt) const override;
        //double eta_v(double b, Vector3d &n_los, double nu, double n_nt) const override;

    protected:
        double s_;
        double gamma_min_;
    private:
        std::string plasma_;
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
        double _nf(const Vector3d &point, double psi) const override;

    private:
        double scale_factor_;
        SimulationInterpolater interp_;
};


class ConstrainedSimulationNField: public PowerLawNField {
    public:
        ConstrainedSimulationNField(Delaunay_triangulation *tr_nt,
                                    Delaunay_triangulation *tr_cold,
                                    bool in_plasma_frame, double s, double gamma_min, bool changing_s=false,
                                    double scale_factor_nt=1.0,
                                    double scale_factor_cold=1.0);
        double _nf(const Vector3d &point, double psi) const override;

    private:
        // Scale factor for NT particles density. Can be arbitrary (e.g. in ``j^2`` heating model), can be (0, 1)-bounded
        // (e.g. in ``n`` or ``B^2`` model)
        double scale_factor_nt_;
        // (0, 1) Maximal fraction of the cold particles that can be heated.
        double scale_factor_cold_;
        SimulationInterpolater interp_nt_;
        SimulationInterpolater interp_cold_;
};


class ConstrainedSigmaSimulationNField: public PowerLawNField {
    public:
        ConstrainedSigmaSimulationNField(Delaunay_triangulation *tr_nt,
                                         Delaunay_triangulation *tr_cold,
                                         Delaunay_triangulation *tr_sigma,
                                         bool in_plasma_frame, double s, double gamma_min, bool changing_s=false,
                                         double scale_factor_nt=1.0,
                                         double scale_factor_cold=1.0);
        double _nf(const Vector3d &point, double psi) const override;

    private:
        // Scale factor for NT particles density. Can be arbitrary (e.g. in ``j^2`` heating model), can be (0, 1)-bounded
        // (e.g. in ``n`` or ``B^2`` model)
        double scale_factor_nt_;
        // (0, 1) Maximal fraction of the cold particles that can be heated.
        double scale_factor_cold_;
        SimulationInterpolater interp_nt_;
        SimulationInterpolater interp_cold_;
        SimulationInterpolater interp_sigma_;
};


class ConstrainedBetaSimulationNField: public PowerLawNField {
    public:
        ConstrainedBetaSimulationNField(Delaunay_triangulation *tr_cold,
                                        Delaunay_triangulation *tr_Bsq,
                                        Delaunay_triangulation *tr_jsq,
                                        Delaunay_triangulation *tr_sigma,
                                        std::string  particles_heating_model,
                                        std::string constrain_type,
                                        bool in_plasma_frame, double s, double gamma_min, bool changing_s=false,
                                        double scale_factor_nt=1.0,
                                        double max_frac_cold=1.0,
                                        double scale_factor_border=0.0,
                                        double psi_mean = 1.0,
                                        double psi_width = 0.1);
        double _nf(const Vector3d &point, double psi) const override;

    private:
        std::string particles_heating_model_;
        std::string constrain_type_;
        // Scale factor for any of the model heating mechanisms NT particles density. Can be arbitrary (e.g. in ``jsq``
        // or ``byhand`` heating model), can be (0, 1)-bounded (e.g. in ``n`` or ``bsq`` model)
        double scale_factor_nt_;
        // Scale factor for NT particles density at the border (0, 1). It modifies underlying NT-density from one of the
        // model heating mechanism.
        double scale_factor_border_;
        // (0, 1) Maximal fraction of the cold particles that can be heated.
        double max_frac_cold_;
        // Psi of border NT particles density enhancement: mean (0, 1) and width of Gaussian profile
        double psi_mean_;
        double psi_width_;
        SimulationInterpolater interp_cold_;
        SimulationInterpolater interp_Bsq_;
        SimulationInterpolater interp_jsq_;
        SimulationInterpolater interp_sigma_;

};


class ByHandSimulationNField: public PowerLawNField {
    public:
        ByHandSimulationNField(Delaunay_triangulation *tr_cold,
                               bool in_plasma_frame, double s, double gamma_min, bool changing_s=false,
                               double max_frac_cold=1.0,
                               double scale_factor_border=1.0,
                               double psi_mean = 1.0,
                               double psi_width = 0.1,
                               double constant_floor_scale = 0.0,
                               double n = 1.0);
        double _nf(const Vector3d &point, double psi) const override;

    private:
        // Scale factor for NT particles density at the border (0, 1).
        double scale_factor_border_;
        // (0, 1) Maximal fraction of the cold particles that can be heated.
        double max_frac_cold_;
        // Psi of border NT particles density enhancement: mean (0, 1) and width of Gaussian profile
        double psi_mean_;
        double psi_width_;
        // Constant additive number density amplitude
        double constant_floor_scale_;
        double n_;
        SimulationInterpolater interp_cold_;
};

#endif //MHD_TRANSFER_NFIELD_H
