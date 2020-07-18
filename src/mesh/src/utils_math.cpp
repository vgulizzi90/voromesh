// utils_math.cpp

#include "utils_math.hpp"

extern "C"
{
#include "../triangle/tricall.h"
#include "../triangle/triangle.h"
}

#include "../tetgen/tetgen.h"

namespace voromesh
{
namespace math
{

// DELAUNAY TRIANGULATION #############################################
double triangle_quality(const double * A, const double * B, const double * C)
{
    const double AB[2] = {B[0]-A[0], B[1]-A[1]};
    const double BC[2] = {B[0]-C[0], B[1]-C[1]};
    const double AC[2] = {C[0]-A[0], C[1]-A[1]};
    const double a = std::sqrt(AB[0]*AB[0]+AB[1]*AB[1]);
    const double b = std::sqrt(BC[0]*BC[0]+BC[1]*BC[1]);
    const double c = std::sqrt(AC[0]*AC[0]+AC[1]*AC[1]);
    const double d = a*b*c;
    const double q = (b+c-a)*(c+a-b)*(a+b-c)/d;

    return q;
}

void delaunay_2d(std::vector<double> & points, std::vector<int> & tri, std::vector<int> & tri_edges,
                 const double min_quality)
{
    // PARAMETERS -----------------------------------------------------
    const int n_points = points.size()/2;
    // ----------------------------------------------------------------

    // VARIABLES ------------------------------------------------------
    int n_edges, n_triangles;
    double q;

    // TRIANGLE
    double * points_ptr = points.data();
    int * edges, * edge_markers, * triangles;
    // ----------------------------------------------------------------

    // CALL TRIANGLE --------------------------------------------------
    tricall(n_points, &points_ptr, &n_edges, &edges, &edge_markers, &n_triangles, &triangles);
    // ----------------------------------------------------------------

    // COPY EDGES AND TRIANGLES FOR OUTPUT ----------------------------
    tri.clear();
    for (int t = 0; t < n_triangles; ++t)
    {
        const double * A = &points[2*triangles[3*t+0]];
        const double * B = &points[2*triangles[3*t+1]];
        const double * C = &points[2*triangles[3*t+2]];

        q = triangle_quality(A, B, C);
        
        if (q > min_quality)
        {
            tri.push_back(triangles[3*t+0]);
            tri.push_back(triangles[3*t+1]);
            tri.push_back(triangles[3*t+2]);
        }
    }

    tri_edges.clear();
    for (int ed = 0; ed < n_edges; ++ed)
    {
        tri_edges.push_back(edges[2*ed+0]);
        tri_edges.push_back(edges[2*ed+1]);
    }
    // ----------------------------------------------------------------

    // FREE MEMORY ----------------------------------------------------
    delete [] edges;
    delete [] edge_markers;
    delete [] triangles;
    // ----------------------------------------------------------------

}

void delaunay_3d(std::vector<double> & points, std::vector<int> & tet, std::vector<int> & tet_edges,
                 const double min_quality)
{
    // PARAMETERS -----------------------------------------------------
    const int n_points = points.size()/3;
    // ----------------------------------------------------------------

    // VARIABLES ------------------------------------------------------
    int n_edges, n_tetrahedra;

    // TETGEN
    tetgenio in, out;
    // ----------------------------------------------------------------

    // INITIALIZATION -------------------------------------------------
    in.numberofpoints = n_points;

    in.pointlist = new double[3*in.numberofpoints];
    std::copy(points.begin(), points.end(), in.pointlist);
    // ----------------------------------------------------------------

    // CALL TETGEN ----------------------------------------------------
    char sw[] = "Qe";

    tetrahedralize(sw, &in, &out);
    // ----------------------------------------------------------------

    // COPY EDGES AND TETRAHEDRA FOR OUTPUT ---------------------------
    tet_edges.resize(2*out.numberofedges);
    tet.resize(4*out.numberoftetrahedra);

    std::copy(out.edgelist, out.edgelist+2*out.numberofedges, tet_edges.data());
    std::copy(out.tetrahedronlist, out.tetrahedronlist+4*out.numberoftetrahedra, tet.data());
    // ----------------------------------------------------------------
}
// ####################################################################

// MATRIX OPERATION ###################################################
namespace matrix
{

void print_inline(const int N, const double * S)
{
    for (int n = 0; n < (N-1); ++n)
    {
        std::cout << std::scientific << std::setprecision(5) << std::setw(12) << S[n] << " ";
    }
    std::cout << std::scientific << std::setprecision(5) << std::setw(12) << S[N-1];
}

void print_inline_int(const int N, const int * S)
{
    for (int n = 0; n < (N-1); ++n)
    {
        std::cout << S[n] << " ";
    }
    std::cout << S[N-1];
}
// ####################################################################
}

}
}