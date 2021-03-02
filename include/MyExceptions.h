#ifndef MHD_TRANSFER_INCLUDE_MYEXCEPTIONS_H
#define MHD_TRANSFER_INCLUDE_MYEXCEPTIONS_H

#include <exception>


class AbsentDataFile : public std::exception {
    std::string msg;
    public:
        AbsentDataFile(const std::string& fname) : msg(std::string("Absent data file : ") + fname) {}
        const char * what () const noexcept override {
            return msg.c_str();
        }
};


class NotImplmentedParticlesHeating : public std::exception {
        std::string msg;
    public:
        NotImplmentedParticlesHeating(const std::string& model) : msg(std::string("Not implemented particles heating model: ") + model) {}
        const char * what () const noexcept override {
            return msg.c_str();
        }
};


class PhysicalException : public std::exception {
        std::string msg;
    public:
        PhysicalException(const std::string& quantity) : msg(std::string("Problems with quantity : ") + quantity) {}
        const char * what () const noexcept override {
            return msg.c_str();
        }
};

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
