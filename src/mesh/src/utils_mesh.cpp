// utils_mesh.cpp

#include "utils_io.hpp"
#include "utils_math.hpp"
#include "utils_mesh.hpp"

#include "../tetgen/tetgen.h"

namespace voromesh
{
namespace utils
{
// DISTANCE OF POINT FROM A SEGMENT ###################################
double d_point_segment(const double * A, const double * B, const double * P)
{
    // PARAMETERS
    const double AB[2] = {B[0]-A[0], B[1]-A[1]};
    const double AP[2] = {P[0]-A[0], P[1]-A[1]};
    const double AP2 = AB[0]*AB[0]+AB[1]*AB[1];
    const double AB_AP = AB[0]*AP[0]+AB[1]*AP[1];
    double aux[2];

    if (AB_AP <= 0.0)
    {
        return std::sqrt(AP2);
    }
    else if (AB_AP >= AP2)
    {
        aux[0] = P[0]-B[0];
        aux[1] = P[1]-B[1];
        return std::sqrt(aux[0]*aux[0]+aux[1]*aux[1]);
    }
    else
    {
        aux[0] = P[0]-(A[0]+AB[0]*AB_AP/AP2);
        aux[1] = P[1]-(A[1]+AB[1]*AB_AP/AP2);
        return std::sqrt(aux[0]*aux[0]+aux[1]*aux[1]);
    }
}
// ####################################################################

// SIGNED DISTANCE OF A CLOUD OF 2D POINTS FROM A POLYGON #############
void sd_points_polygon(const std::vector<double> & vertices,
                       const std::vector<double> & points,
                       std::vector<double> & d)
{
    // PARAMETERS
    const int n_vertices = vertices.size()/2;
    const int n_points = points.size()/2;

    // RESIZE MEMORY
    d.resize(n_points, std::numeric_limits<double>::max());

    // LOOP OVER THE POINTS
    for (int p = 0; p < n_points; ++p)
    {
        const double * P = &points[2*p];

        int in = 0;

        for (int v = 0; v < n_vertices; ++v)
        {
            const double * Vi = ((v == 0) ? &vertices[2*(n_vertices-1)] : &vertices[2*(v-1)]);
            const double * Vj = &vertices[2*v];

            if (Vi[1] <= P[1])
            {
                if ((Vj[1] > P[1]) && ((((Vj[0]-Vi[0])*(P[1]-Vi[1])-(Vj[1]-Vi[1])*(P[0]-Vi[0]))) > 0.0))
                {
                    in += 1;
                }
            }
            else
            {
                if ((Vj[1] <= P[1]) && ((((Vj[0]-Vi[0])*(P[1]-Vi[1])-(Vj[1]-Vi[1])*(P[0]-Vi[0]))) < 0.0))
                {
                    in -= 1;
                }
            }

            d[p] = std::min(d[p], d_point_segment(Vi, Vj, P));
        }

        if (in != 0)
        {
            d[p] *= -1.0;
        }
    }
}
// ####################################################################


// COMPUTE THE MESH OF A POLYGON ######################################
void mesh_polygon(const double mesh_size,
                  const std::vector<double> & vertices,
                  const std::vector<double> & boundary_nodes,
                  std::vector<double> & nodes, std::vector<int> & triangles)
{
    // PARAMETERS -----------------------------------------------------
    const int n_vertices = vertices.size()/2;
    const int n_boundary_nodes = boundary_nodes.size()/2;

    // DISTMESH ALGORITHM
    const double dX_tol = 0.001;
    const double ttol = 0.1;
    const double Fscale = 1.2;
    const double dt = 0.2;
    const double geps = 0.001*mesh_size;
    const double deps = std::sqrt(std::numeric_limits<double>::epsilon())*mesh_size;
    const double big_eps = mesh_size*std::sqrt(3.0)/32.0;
    
    const int max_it = 100;
    const int control_frequency = 30;
    // ----------------------------------------------------------------

    // VARIABLES ------------------------------------------------------
    double X1_ends[2], X2_ends[2];
    std::vector<double> dpp, tri_c;
    int n_nodes, n_triangles, n_edges, n_short_edges;
    std::vector<int> edges, short_edges, nodes_to_be_removed;
    std::vector<int>::iterator iter;
    std::vector<double> edges_vec, edges_h;
    std::vector<double> Fn, dX;
    double dX_max;
    // ----------------------------------------------------------------

    // BOUNDING BOX ---------------------------------------------------
    X1_ends[0] = vertices[2*0+0]; X1_ends[1] = vertices[2*0+0];
    X2_ends[0] = vertices[2*0+1]; X2_ends[1] = vertices[2*0+1];
    for (int v = 1; v < n_vertices; ++v)
    {
        X1_ends[0] = std::min(X1_ends[0], vertices[2*v+0]);
        X1_ends[1] = std::max(X1_ends[1], vertices[2*v+0]);
        X2_ends[0] = std::min(X2_ends[0], vertices[2*v+1]);
        X2_ends[1] = std::max(X2_ends[1], vertices[2*v+1]);
    }

    /* DEBUG
    std::cout << "vertices:" << std::endl;
    for (int n = 0; n < n_vertices; ++n)
    {
        math::matrix::print_inline(2, &vertices[2*n]); std::cout << std::endl;
    }
    std::cout << "bounding box:" << std::endl;
    std::cout << "X1_ends: " << X1_ends[0] << "," << X1_ends[1] << std::endl;
    std::cout << "X2_ends: " << X2_ends[0] << "," << X2_ends[1] << std::endl;
    */
    // ----------------------------------------------------------------

    // COPY THE BOUNDARY NODES ----------------------------------------
    nodes.resize(2*n_boundary_nodes);
    std::copy(boundary_nodes.begin(), boundary_nodes.end(), nodes.begin());
    // ----------------------------------------------------------------

    // ADD AN INITIAL GRID-LIKE DISTRIBUTION OF NODES -----------------
    {
        const int n1 = round((X1_ends[1]-X1_ends[0])/mesh_size);
        const int n2 = round((X2_ends[1]-X2_ends[0])/(mesh_size*0.5*std::sqrt(3.0)));

        if ((n1 > 1) && (n2 > 1))
        {
            std::vector<double> X1(n1, 0.0);
            std::vector<double> X2(n2, 0.0);
            math::linspace(X1_ends[0], X1_ends[1], n1, X1.data());
            math::linspace(X2_ends[0], X2_ends[1], n2, X2.data());

            const double dX1 = (X1_ends[1]-X1_ends[0])/(n1-1);

            for (int i = 0; i < n1; ++i)
            for (int j = 0; j < n2; ++j)
            {
                nodes.push_back(X1[i]+(j%2)*0.5*dX1);
                nodes.push_back(X2[j]);
            }
        }
    }
    n_nodes = nodes.size()/2;

    /* DEBUG
    std::cout << "nodes:" << std::endl;
    for (int n = 0; n < n_nodes; ++n)
    {
        math::matrix::print_inline(2, &nodes[2*n]); std::cout << std::endl;
    }
    */
    // ----------------------------------------------------------------

    // REMOVE THE NODES OUTSIDE THE POLYGON ---------------------------
    sd_points_polygon(vertices, nodes, dpp);

    for (int n = (n_nodes-1); n >= n_boundary_nodes; --n)
    {
        if (dpp[n] > -geps)
        {
            nodes.erase(nodes.begin()+2*n+1);
            nodes.erase(nodes.begin()+2*n);
        }
    }
    n_nodes = nodes.size()/2;

    /* DEBUG
    std::cout << "distance:" << std::endl;
    for (int n = 0; n < n_nodes; ++n)
    {
        std::cout << dpp[n] << std::endl;
    }
    */

    /* DEBUG
    std::cout << "nodes:" << std::endl;
    for (int n = 0; n < n_nodes; ++n)
    {
        math::matrix::print_inline(2, &nodes[2*n]); std::cout << std::endl;
    }
    */
    // ----------------------------------------------------------------

    // START THE ITERATIVE ALGORITHM ----------------------------------
    int it = 0;
    bool flag_continue = true;

    while ((it < max_it) and flag_continue)
    {
        // TRIANGULATE NODES
        math::delaunay_2d(nodes, triangles, edges);
        n_triangles = triangles.size()/3;
        n_edges = edges.size()/2;

        /* DEBUG
        std::cout << "it: " << it << ",nodes:" << std::endl;
        for (int n = 0; n < n_boundary_nodes; ++n)
        {
            math::matrix::print_inline(2, &nodes[2*n]); std::cout << std::endl;
        }
        std::cout << std::endl;
        for (int n = n_boundary_nodes; n < n_nodes; ++n)
        {
            math::matrix::print_inline(2, &nodes[2*n]); std::cout << std::endl;
        }
        std::cout << "triangles:" << std::endl;
        for (int n = 0; n < triangles.size()/3; ++n)
        {
            math::matrix::print_inline_int(3, &triangles[3*n]); std::cout << std::endl;
        }
        */

        // REMOVE TRIANGLES THAT FALL OUTSITE
        tri_c.resize(2*n_triangles);
        for (int t = 0; t < n_triangles; ++t)
        {
            math::tri_centroid_2d(&nodes[2*triangles[3*t+0]],
                                  &nodes[2*triangles[3*t+1]],
                                  &nodes[2*triangles[3*t+2]],
                                  &tri_c[2*t]);
        }
        sd_points_polygon(vertices, tri_c, dpp);
        for (int t = (n_triangles-1); t >= 0; --t)
        {
            if (dpp[t] > geps)
            {
                triangles.erase(triangles.begin()+3*t+2);
                triangles.erase(triangles.begin()+3*t+1);
                triangles.erase(triangles.begin()+3*t);
            }
        }
        n_triangles = triangles.size()/3;

        // EDGE LENGTHS
        edges_vec.resize(2*n_edges);
        edges_h.resize(n_edges);
        for (int ed = 0; ed < n_edges; ++ed)
        {
            edges_vec[2*ed+0] = nodes[2*edges[2*ed+1]+0]-nodes[2*edges[2*ed]+0];
            edges_vec[2*ed+1] = nodes[2*edges[2*ed+1]+1]-nodes[2*edges[2*ed]+1];
            edges_h[ed] = std::sqrt(edges_vec[2*ed+0]*edges_vec[2*ed+0]+edges_vec[2*ed+1]*edges_vec[2*ed+1]);
        }

        // LOOK FOR SHORT EDGES
        short_edges.clear();
        for (int ed = 0; ed < n_edges; ++ed)
        {
            if (edges_h[ed] < 0.5*mesh_size)
            {
                short_edges.push_back(ed);
            }
        }
        n_short_edges = short_edges.size();

        if (control_frequency%(it+1) == 0)
        {
            nodes_to_be_removed.clear();
            for (int ed = 0; ed < n_short_edges; ++ed)
            {
                if (edges[2*short_edges[ed]+0] >= n_boundary_nodes)
                {
                    nodes_to_be_removed.push_back(edges[2*short_edges[ed]+0]);
                }
                if (edges[2*short_edges[ed]+1] >= n_boundary_nodes)
                {
                    nodes_to_be_removed.push_back(edges[2*short_edges[ed]+1]);
                }
            }
            std::sort(nodes_to_be_removed.begin(), nodes_to_be_removed.end());
            iter = std::unique(nodes_to_be_removed.begin(), nodes_to_be_removed.end());
            nodes_to_be_removed.erase(iter, nodes_to_be_removed.end());

            for (int i = (nodes_to_be_removed.size()-1); i >= 0; --i)
            {
                int n = nodes_to_be_removed[i];
                nodes.erase(nodes.begin()+2*n+1);
                nodes.erase(nodes.begin()+2*n);
            }
            n_nodes = nodes.size()/2;

            it += 1;
            continue;
        }

        // COMPUTE THE NODAL FORCES
        Fn.resize(2*n_nodes, 0.0);
        for (int ed = 0; ed < n_edges; ++ed)
        {
            const double F_edge_0 = (std::max(mesh_size-edges_h[ed], 0.0)/edges_h[ed])*edges_vec[2*ed+0];
            const double F_edge_1 = (std::max(mesh_size-edges_h[ed], 0.0)/edges_h[ed])*edges_vec[2*ed+1];

            Fn[2*edges[2*ed+0]+0] -= F_edge_0;
            Fn[2*edges[2*ed+0]+1] -= F_edge_1;
            Fn[2*edges[2*ed+1]+0] += F_edge_0;
            Fn[2*edges[2*ed+1]+1] += F_edge_1;
        }
        // Boundary nodes shall not move
        for (int n = 0; n < n_boundary_nodes; ++n)
        {
            Fn[2*n+0] = 0.0;
            Fn[2*n+1] = 0.0;
        }

        // MOVE NODES
        dX.resize(2*n_nodes, 0.0);
        for (int n = n_boundary_nodes; n < n_nodes; ++n)
        {
            dX[2*n+0] = dt*Fn[2*n+0];
            dX[2*n+1] = dt*Fn[2*n+1];
            nodes[2*n+0] += dX[2*n+0];
            nodes[2*n+1] += dX[2*n+1];
        }

        // MAKE SURE NO NODES OTHER THAN THE BOUNDARY NODES LIE ON THE
        // BOUNDARY
        /*
        sd_points_polygon(vertices, nodes, dpp);
        for (int n = (n_nodes-1); n >= n_boundary_nodes; --n)
        {
            if (dpp[n] > (-big_eps))
            {
                nodes.erase(nodes.begin()+2*n+1);
                nodes.erase(nodes.begin()+2*n);
            }
        }
        n_nodes = nodes.size()/2;
        */

        // CHECK THE MAXIMUM DISPLACEMENT
        dX_max = 0.0;
        for (int n = n_boundary_nodes; n < n_nodes; ++n)
        {
            dX_max = std::max(dX_max, std::sqrt(dX[2*n+0]*dX[2*n+0]+dX[2*n+1]*dX[2*n+1]));
        }
        flag_continue = (dX_max > dX_tol);

        // NEXT ITERATION
        it += 1;
    }

    math::delaunay_2d(nodes, triangles, edges);
    // ----------------------------------------------------------------

    /* DEBUG
    std::cout << "nodes:" << std::endl;
    for (int n = 0; n < n_nodes; ++n)
    {
        math::matrix::print_inline(2, &nodes[2*n]); std::cout << std::endl;
    }
    */
    /* DEBUG
    std::cout << "tri_c:" << std::endl;
    for (int n = 0; n < n_triangles; ++n)
    {
        math::matrix::print_inline(2, &tri_c[2*n]); std::cout << std::endl;
    }
    */
    /* DEBUG
    std::cout << "triangles:" << std::endl;
    for (int n = 0; n < triangles.size()/3; ++n)
    {
        math::matrix::print_inline_int(3, &triangles[3*n]); std::cout << std::endl;
    }
    */
}
// ####################################################################

// COMPUTE THE MESH OF A POLYHEDRON ###################################
void c_mesh_convex_polyhedron(const double * centroid,
                              const std::vector<double> & boundary_nodes, const std::vector<int> & boundary_triangles,
                              std::vector<double> & nodes, std::vector<int> & tetrahedra)
{
    // PARAMETERS -----------------------------------------------------
    const int n_boundary_nodes = boundary_nodes.size()/3;
    const int n_boundary_triangles = boundary_triangles.size()/3;
    // ----------------------------------------------------------------

    // MAKE THE CENTROIDAL MESH ---------------------------------------
    // We construct the tetrahedra by connecting each boundary triangle
    // with the centroid of the polyhedron.
    // ----------------------------------------------------------------
    nodes.resize(3*(n_boundary_nodes+1));
    tetrahedra.resize(4*n_boundary_triangles);

    // NODES
    std::copy(boundary_nodes.begin(), boundary_nodes.end(), nodes.begin());
    nodes[3*n_boundary_nodes+0] = centroid[0];
    nodes[3*n_boundary_nodes+1] = centroid[1];
    nodes[3*n_boundary_nodes+2] = centroid[2];

    // TETRAHEDRA
    for (int t = 0; t < n_boundary_triangles; ++t)
    {
        tetrahedra[4*t+0] = boundary_triangles[3*t+0];
        tetrahedra[4*t+1] = boundary_triangles[3*t+1];
        tetrahedra[4*t+2] = boundary_triangles[3*t+2];
        tetrahedra[4*t+3] = n_boundary_nodes;
    }
    // ----------------------------------------------------------------
}

void mesh_convex_polyhedron(const double mesh_size, const std::vector<std::array<double, 4>> & planes,
                            const std::vector<double> & boundary_nodes, const std::vector<int> & boundary_triangles,
                            std::vector<double> & nodes, std::vector<int> & tetrahedra)
{
    // PARAMETERS -----------------------------------------------------
    const int n_boundary_nodes = boundary_nodes.size()/3;
    // ----------------------------------------------------------------

    // VARIABLES ------------------------------------------------------
    // TETGEN
    tetgenio in, out;
    // ----------------------------------------------------------------

    // INITIALIZE THE TETGEN INPUT ------------------------------------
    // POINTS
    in.numberofpoints = n_boundary_nodes;
    in.pointlist = new double[3*in.numberofpoints];
    std::copy(boundary_nodes.begin(), boundary_nodes.end(), in.pointlist);

    // FACES
    // ----------------------------------------------------------------
}
// ####################################################################


}
}