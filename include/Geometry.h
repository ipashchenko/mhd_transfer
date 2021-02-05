#ifndef MHD_TRANSFER_GEOMETRY_H
#define MHD_TRANSFER_GEOMETRY_H

#include <Eigen/Eigen>
#include "Ray.h"
#include "Intersection.h"
#include <list>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <list>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

using Eigen::Vector3d;

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef K::Ray_3 Ray_3;
typedef K::Line_3 Line_3;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Primitive_id Primitive_id;


std::list<double> intersection(Vector3d R0, Vector3d Rd, double A = 0, double B = 0, double C = 0, double D = 0,
                               double E = 0, double F = 0, double G = 0, double H = 0, double I = 0, double J = 0);

class Ray;
class Intersection;

class BaseGeometry {
    public:
        virtual std::list<Intersection> hit(Ray &ray) const = 0;
        virtual bool is_within(Vector3d& point) const = 0;
};

class Geometry : public BaseGeometry {
    public:
        virtual const double big_scale() const = 0;
        virtual const Vector3d& origin() const = 0;
        virtual double radius_at_given_distance(const Vector3d& point) const = 0;
        std::pair<Vector3d, Vector3d> full_infinite_path(Ray &ray) const;
        std::pair<Vector3d, Vector3d> half_infinite_path(Ray &ray, const Vector3d &point) const;

};


// TODO: Simulation geometry could be defined by analytical y = f(x) dependence of radius on the distance z!
class SimulationGeometry : public BaseGeometry {
    public:
        explicit SimulationGeometry(Tree *tree);
        std::list<Intersection> hit(Ray &ray) const override;
        double radius_at_given_distance(const Vector3d& point) const;
        bool is_within(Vector3d& point) const;

    private:
        Tree* tree_;
};


#endif //MHD_TRANSFER_GEOMETRY_H
