#ifndef MHD_TRANSFER_SIMULATIONINTERPOLATER_H
#define MHD_TRANSFER_SIMULATIONINTERPOLATER_H


#include <Eigen/Eigen>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>

typedef CGAL::Simple_cartesian<double>                                   K_;
typedef K_::Point_2                                                Point_;
typedef CGAL::Triangulation_vertex_base_with_info_2<double, K_>      Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                  Tds;
typedef CGAL::Delaunay_triangulation_2<K_, Tds>                    Delaunay_triangulation;
typedef K_::FT                                               Coord_type;
typedef std::vector<Coord_type >                            Scalar_vector;
typedef CGAL::Barycentric_coordinates::Triangle_coordinates_2<K_> Triangle_coordinates;

using Eigen::Vector2d;

// Using file with (r, z, value) or (r, Psi, value) create triangulation. r & z are in pc
void create_triangulation(std::string fn, Delaunay_triangulation* tr);


class SimulationInterpolater {
    public:
        explicit SimulationInterpolater(Delaunay_triangulation *tr, double nan_value=0);
        // Interpolate in (r, z) or (Psi, z) plane, where r & z are in pc
        double interpolated_value(Vector2d point) const;
        // Find gradient in (r, z) coordinates only!
        Vector2d gradient(Vector2d point) const;

    private:
        double nan_value_;
        Delaunay_triangulation* tr_;
        mutable Delaunay_triangulation::Face_handle previous_hit_fh_;
};


#endif //MHD_TRANSFER_SIMULATIONINTERPOLATER_H
