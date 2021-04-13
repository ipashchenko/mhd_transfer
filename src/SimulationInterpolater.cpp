#include "SimulationInterpolater.h"
#include "utils.h"

#define CGAL_HAS_THREADS


SimulationInterpolater::SimulationInterpolater(Delaunay_triangulation *tr, double nan_value) {
    tr_ = tr;
    std::cout << "Initializing interpolater with nan_value = " << nan_value << std::endl;
    nan_value_ = nan_value;
    previous_hit_fh_ = nullptr;

//    for (Delaunay_triangulation::Finite_faces_iterator fit = tr->finite_faces_begin(); fit != tr->finite_faces_end(); ++fit) {
//        Delaunay_triangulation::Face_handle fh = fit;
//        previous_hit_fh_ = std::make_shared<Delaunay_triangulation::Face_handle>(fh);
//        break;
//    }
}

double SimulationInterpolater::interpolated_value(Vector2d point) const {
    Point_ pt(point[0], point[1]);

    Delaunay_triangulation::Face_handle fh;
    // TODO: Try T.inexact_locate for speed
    if (previous_hit_fh_ != nullptr) {
//        std::cout << "Hit with hint" << std::endl;
        fh = tr_->locate(pt, previous_hit_fh_);
//        std::cout << "Done Hit with hint" << std::endl;
    } else {
//        std::cout << "First time hit" << std::endl;
        fh = tr_->locate(pt);
    }

    // This previous w/o hints
//    fh = tr_->locate(pt);

    if (tr_->is_infinite(fh)) {
        previous_hit_fh_ = nullptr;
        return nan_value_;
    } else {
        previous_hit_fh_ = fh;
    }

    std::vector<Point_ > vertexes;
    std::vector<double> info;

    for (int i=0; i<3; i++) {
        vertexes.push_back(fh->vertex(i)->point());
        info.push_back(fh->vertex(i)->info());

//        std::cout << "Triangle:\t" << tr_->triangle(fh) << std::endl;
//        std::cout << "Vertex 0:\t" << tr_->triangle(fh)[i] << std::endl;
//        std::cout << "Value:\t" << fh->vertex(i)->info() << std::endl;
    }

    // Create an std::vector to store coordinates.
    Scalar_vector coordinates;
    // Instantiate class Triangle_coordinates_2 for the triangle defined above.
    Triangle_coordinates triangle_coordinates(vertexes[0], vertexes[1], vertexes[2]);
    triangle_coordinates(pt, std::inserter(coordinates, coordinates.end()));

    double interpolated_value = 0;
    for(int j = 0; j < 3; ++j) {
        interpolated_value += coordinates[j]*info[j];
    }

    if (std::isnan(interpolated_value)) {
        interpolated_value = nan_value_;
    }

    return interpolated_value;
}


Vector2d SimulationInterpolater::gradient(Vector2d point, const SimulationInterpolater& psi_space_interpolater) const {
    double eps_r = 0.025;
    double eps_z = 0.1;
    // point = (r_p, z)
    double r_p = point[0];
    double z = point[1];

//    std::cout << "Calculating gradient in point (r_p, z) = " << r_p << ", " << z << "\n";

    // d/dr
    double d_dr, psi_plus, psi_minus, psi_cur;
//    double h = (1. + r_p)*sqrt(eps);
    double h = r_p*eps_r;
//    std::cout << "h for d/dr = " << h << "\n";
    if(r_p - h <= 0) {
        // Find (Psi, z) coordinate for (r_p+h, z)
        psi_plus = interpolated_value({r_p + h, z});
        psi_cur = interpolated_value({r_p, z});
        d_dr = (psi_space_interpolater.interpolated_value({psi_plus, z}) -
                psi_space_interpolater.interpolated_value({psi_cur, z}))/h;
//        std::cout << "Border d/dr = " << d_dr << "\n";
    } else {
        // Find (Psi, z) coordinates for (r_p+h, z) and (r_p-h, z)
        psi_plus = interpolated_value({r_p + h, z});
        psi_minus = interpolated_value({r_p - h, z});

        d_dr = (psi_space_interpolater.interpolated_value({psi_plus, z}) -
                psi_space_interpolater.interpolated_value({psi_minus, z}))/(2.0*h);
//        std::cout << "Two sided d/dr = " << d_dr << "\n";
    }

    // d/dz
    double d_dz;
//    h = (1. + z)*sqrt(eps);
    h = z*eps_z;
//    std::cout << "h for d/dz = " << h << "\n";
    // Find (Psi, z) coordinates for (r_p, z+h) and (r_p, z-h)
    psi_plus = interpolated_value({r_p, z + h});
    psi_minus = interpolated_value({r_p, z - h});
    if(psi_minus == 1.0){
        psi_cur = interpolated_value({r_p, z});
        d_dz = (psi_space_interpolater.interpolated_value({psi_plus, z + h}) -
                psi_space_interpolater.interpolated_value({psi_cur, z}))/h;
    } else {
        d_dz = (psi_space_interpolater.interpolated_value({psi_plus, z + h}) -
                psi_space_interpolater.interpolated_value({psi_minus, z - h}))/(2.0*h);
    }

//    std::cout << "d/dz = " << d_dz << "\n";
    return {d_dr, d_dz};
}


void create_triangulation(std::string fn, Delaunay_triangulation *tr) {
    std::vector< std::vector<double> > all_points;
    read_from_txt(fn, all_points);
    size_t nrows = all_points.size();

    std::vector< std::pair<Point_,double> > points;
    for (int i=0; i<nrows; i++) {
//        Point_ pt(all_points[i][0], all_points[i][1]);
        // (r, z) or (Psi, z)
        Point_ pt(all_points[i][1], all_points[i][0]);
        points.emplace_back( std::make_pair( pt,  all_points[i][2]) );
    }
    tr->insert(points.begin(), points.end());
}
