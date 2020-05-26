// utils_mesh.hpp

#ifndef VOROMESH_UTILS_HPP_
#define VOROMESH_UTILS_HPP_

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

namespace voromesh
{
namespace utils
{
// DISTANCE OF POINT FROM A SEGMENT ###################################
double d_point_segment(const double * A, const double * B, const double * P);
// ####################################################################


// SIGNED DISTANCE OF A CLOUD OF 2D POINTS FROM A POLYGON #############
void sd_points_polygon(const std::vector<double> & vertices,
                       const std::vector<double> & points,
                       std::vector<double> & d);
// ####################################################################

// COMPUTE THE MESH OF A POLYGON ######################################
void mesh_polygon(const double mesh_size,
                  const std::vector<double> & vertices,
                  const std::vector<double> & boundary_nodes,
                  std::vector<double> & nodes, std::vector<int> & triangles);
// ####################################################################

// COMPUTE THE MESH OF A POLYHEDRON ###################################
void c_mesh_convex_polyhedron(const double * centroid,
                              const std::vector<double> & boundary_nodes, const std::vector<int> & boundary_triangles,
                              std::vector<double> & nodes, std::vector<int> & tetrahedra);

void mesh_convex_polyhedron(const double mesh_size, const std::vector<std::array<double, 4>> & planes,
                            const std::vector<double> & boundary_nodes, const std::vector<int> & boundary_triangles,
                            std::vector<double> & nodes, std::vector<int> & tetrahedra);
// ####################################################################

}
}

#endif