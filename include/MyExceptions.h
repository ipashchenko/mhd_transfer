#ifndef MHD_TRANSFER_INCLUDE_MYEXCEPTIONS_H
#define MHD_TRANSFER_INCLUDE_MYEXCEPTIONS_H

#include <exception>


class BadParticleProfileParameters : public std::exception {
    const char * what () const noexcept override {
        return "Check particle profile parameters!";
    }
};

class BadPlasmaContentParameters : public std::exception {
    const char * what () const noexcept override {
        return "Check plasma content parameters!";
    }
};

class PointOutsideOfPolyhedron : public std::exception {
    const char * what () const noexcept override {
        return "Point has z-coordinate that lies outside of Simulation Geometry Polyhedron!";
    }
};

#endif //MHD_TRANSFER_INCLUDE_MYEXCEPTIONS_H
