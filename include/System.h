#ifndef MHD_TRANSFER_SYSTEMS_H
#define MHD_TRANSFER_SYSTEMS_H

#include <iostream>
#include "Jet.h"


typedef std::vector<double> state_type;

class System {
    public:
        System(Jet* newjet, Vector3d &newpoint_in, Vector3d &newray_direction, double newnu);
        virtual void operator() (const double &x, double &dxdt, double t) = 0;

	protected:
        Jet* jet;
        Vector3d point_in;
        Vector3d ray_direction;
        double nu;
};

class Tau : public System {
    public:
        Tau(Jet* newjet, Vector3d &newpoint_in, Vector3d &newray_direction, double newnu) :
            System(newjet, newpoint_in,newray_direction, newnu) {};
        void operator() (const double &x, double &dxdt, double t) override;
};


class TauFR : public System {
    public:
        TauFR(Jet* newjet, Vector3d &newpoint_in, Vector3d &newray_direction, double newnu) :
            System(newjet, newpoint_in,newray_direction, newnu) {};
        void operator() (const double &x, double &dxdt, double t) override;
};


// point_in is actually furthest point (i.e. ``point_out`` in intersection).
class I : public System {
 public:
  I(Jet* newjet, Vector3d &newpoint_in, Vector3d &newray_direction,
    double newnu) : System(newjet, newpoint_in, newray_direction, newnu) {};

  void operator() (const double &x, double &dxdt, double t) override;
};


class Speed : public System {
    public:
        Speed(Jet* newjet, Vector3d &newpoint_in, Vector3d &newray_direction,
          double newnu) : System(newjet, newpoint_in, newray_direction, newnu) {};

        void operator() (const double &x, double &dxdt, double t) override;
};

class FullStokes {
    public:
		FullStokes(Jet* newjet, Vector3d &newpoint_in, Vector3d &newray_direction, double newnu);
		void operator() (const state_type &x, state_type &dxdt, double t);

    protected:
		Jet* jet;
		Vector3d point_in;
		Vector3d ray_direction;
		double nu;
};

bool check_opt_depth(double tau_max, const double &x);


#endif //MHD_TRANSFER_SYSTEMS_H
