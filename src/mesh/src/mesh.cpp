// mesh.cpp

#include <limits>
#include "mesh.hpp"
#include "utils_io.hpp"
#include "utils_vtk.hpp"
#include "utils_mesh.hpp"

namespace voromesh
{
// VORONOI VERTEX CLASS ###############################################
// CONSTRUCTOR ========================================================
VoronoiVertex::VoronoiVertex()
:
state(0)
{
    this->x[0] = 0.0;
    this->x[1] = 0.0;
    this->x[2] = 0.0;
}
// ====================================================================
// ####################################################################



// VORONOI EDGE CLASS #################################################
// CONSTRUCTOR ========================================================
VoronoiEdge::VoronoiEdge()
:
state(0)
{

}
// ====================================================================

// BUILD MESH =========================================================
void VoronoiEdge::build_mesh(const double mesh_size,
                             const std::vector<VoronoiVertex> & vertices,
                             std::vector<double> & nodes, std::vector<int> & conn) const
{
    // PARAMETERS -----------------------------------------------------
    const VoronoiVertex & verA = vertices[this->v_conn[0]];
    const VoronoiVertex & verB = vertices[this->v_conn[1]];
    const double * A = verA.x;
    const double * B = verB.x;
    const double edge_len = math::L2_distance(A, B, 3);
    // ----------------------------------------------------------------

    // VARIABLES ------------------------------------------------------
    int n_elems;
    double dV[3];
    // ----------------------------------------------------------------

    if (edge_len > 1.3*mesh_size)
    {
        // NUMBER OF ELEMENTS
        n_elems = ceil(edge_len/mesh_size);
        dV[0] = (B[0]-A[0])/n_elems;
        dV[1] = (B[1]-A[1])/n_elems;
        dV[2] = (B[2]-A[2])/n_elems;

        // RESIZE NODES
        nodes.resize(3*(n_elems+1));

        // RESIZE CONNECTIVITY
        conn.resize(2*n_elems);

        // FIRST ELEMENT
        nodes[3*0+0] = A[0];
        nodes[3*0+1] = A[1];
        nodes[3*0+2] = A[2];
        nodes[3*2+0] = A[0]+dV[0];
        nodes[3*2+1] = A[1]+dV[1];
        nodes[3*2+2] = A[2]+dV[2];

        conn[2*0+0] = 0;
        conn[2*0+1] = 2;

        // ELEMENTS IN THE MIDDLE
        for (int e = 1; e < (n_elems-1); ++e)
        {
            nodes[3*(e+2)+0] = A[0]+(e+1)*dV[0];
            nodes[3*(e+2)+1] = A[1]+(e+1)*dV[1];
            nodes[3*(e+2)+2] = A[2]+(e+1)*dV[2];

            conn[2*e+0] = e+1;
            conn[2*e+1] = e+2;
        }

        // LAST ELEMENT
        nodes[3*1+0] = B[0];
        nodes[3*1+1] = B[1];
        nodes[3*1+2] = B[2];

        conn[2*(n_elems-1)+0] = n_elems;
        conn[2*(n_elems-1)+1] = 1;

        /* DEBUG
        if ((std::abs(dV[1]) < 1.0e-5) and (std::abs(dV[2]) < 1.0e-5))
        {
            std::cout << "dV: "; math::matrix::print_inline(3, dV); std::cout << std::endl;
            std::cout << "A: "; math::matrix::print_inline(3, A); std::cout << std::endl;
            std::cout << "B: "; math::matrix::print_inline(3, B); std::cout << std::endl;
            const int n_nodes = (n_elems+1);
            std::cout << "nodes: " << std::endl;
            for (int n = 0; n < n_nodes; ++n)
            {
                math::matrix::print_inline(3, &nodes[3*n]); std::cout << std::endl;
            }

            io::error("mesh.cpp - VoronoiEdge::build_mesh", "HERE WE ARE");
        }
        */
    }
    else
    {
        // NUMBER OF ELEMENTS
        n_elems = 1;

        // RESIZE NODES
        nodes.resize(3*(n_elems+1));

        // RESIZE CONNECTIVITY
        conn.resize(2*n_elems);

        // FIRST AND LAST ELEMENT
        nodes[3*0+0] = A[0];
        nodes[3*0+1] = A[1];
        nodes[3*0+2] = A[2];
        nodes[3*1+0] = B[0];
        nodes[3*1+1] = B[1];
        nodes[3*1+2] = B[2];

        conn[2*0+0] = 0;
        conn[2*0+1] = 1;
    }
}
// ====================================================================
// ####################################################################



// VORONOI FACE CLASS #################################################
// CONSTRUCTOR ========================================================
VoronoiFace::VoronoiFace()
:
state(0),
bou_type(BOU_TYPE_UNDEFINED)
{
    this->parent_cells[0] = -1;
    this->parent_cells[1] = -1;
}
// ====================================================================

// EVAL LOCAL REFERENCE SYSTEM ========================================
void VoronoiFace::eval_RS(const double * vertices, const int n_vertices, const int * conn)
{
    // PARAMETERS
    const double * V0 = &vertices[3*conn[0]];
    const double * V1 = &vertices[3*conn[1]];

    // VARIABLES
    double aux;
    double M[3];
    double t1[3], t2[3], un[3];

    // EVAL AN AVERAGE POINT AMONG THE FACE VERTICES
    M[0] = 0.0;
    M[1] = 0.0;
    M[2] = 0.0;
    for (int v = 0; v < n_vertices; ++v)
    {
        const double * V = &vertices[3*conn[v]];

        M[0] = (M[0]*v+V[0])/(v+1);
        M[1] = (M[1]*v+V[1])/(v+1);
        M[2] = (M[2]*v+V[2])/(v+1);
    }

    // EVAL THE FIRST UNIT VECTOR
    t1[0] = V1[0]-V0[0];
    t1[1] = V1[1]-V0[1];
    t1[2] = V1[2]-V0[2];
    aux = 1.0/math::L2_norm(t1, 3);
    t1[0] *= aux;
    t1[1] *= aux;
    t1[2] *= aux;

    // EVAL A TEMPORARY SECOND UNIT VECTOR
    t2[0] = M[0]-V0[0];
    t2[1] = M[1]-V0[1];
    t2[2] = M[2]-V0[2];
    aux = 1.0/math::L2_norm(t2, 3);
    t2[0] *= aux;
    t2[1] *= aux;
    t2[2] *= aux;

    // EVAL THE UNIT NORMAL
    math::cross(t1, t2, un);
    aux = 1.0/math::L2_norm(un, 3);
    un[0] *= aux;
    un[1] *= aux;
    un[2] *= aux;

    // EVAL THE SECOND UNIT VECTOR
    math::cross(un, t1, t2);

    // STORE THE UNIT VECTORS
    this->R[0+0*3] = t1[0]; this->R[0+1*3] = t1[1]; this->R[0+2*3] = t1[2];
    this->R[1+0*3] = t2[0]; this->R[1+1*3] = t2[1]; this->R[1+2*3] = t2[2];
    this->R[2+0*3] = un[0]; this->R[2+1*3] = un[1]; this->R[2+2*3] = un[2];
}
// ====================================================================

// BUILD MESH =========================================================
void VoronoiFace::build_mesh(const double mesh_size,
                             const std::vector<VoronoiVertex> & vertices,
                             const std::vector<VoronoiEdge> & edges,
                             const std::vector<double> & global_mesh_nodes,
                             const std::vector<int> & global_mesh_offset,
                             const std::vector<int> & global_mesh_conn,
                             std::vector<double> & nodes, std::vector<int> & conn,
                             std::vector<int> & bou_nodes) const
{
    // PARAMETERS -----------------------------------------------------
    const int n_face_vertices = this->v_conn.size();
    const int n_face_edges = this->ed_conn.size();
    // ----------------------------------------------------------------

    // VARIABLES ------------------------------------------------------
    int e, elem_offset;
    const int * elem_conn;
    int v, ed, n_edge_elems;
    int node_pick;
    int n_boundary_nodes;

    std::vector<double> face_vertices, projected_face_vertices, projected_face_vertices_2d;
    std::vector<double> projected_nodes, projected_nodes_2d;
    std::vector<double> nodes_2d;

    double z, RT[9];
    int n_nodes;
    // ----------------------------------------------------------------

    // CLEAR MEMORY --
    nodes.clear();
    conn.clear();
    bou_nodes.clear();
    // ---------------

    // GATHER THE VERTICES --------------------------------------------
    for (int fv = 0; fv < n_face_vertices; ++fv)
    {
        v = this->v_conn[fv];
        const VoronoiVertex & vertex = vertices[v];

        face_vertices.push_back(vertex.x[0]);
        face_vertices.push_back(vertex.x[1]);
        face_vertices.push_back(vertex.x[2]);
    }

    // PROJECT THE FACE VERTICES ONTO THE FACE LOCAL REF SYSTEM
    projected_face_vertices.resize(face_vertices.size());
    math::matrix::multiply(3, 3, n_face_vertices, this->R, face_vertices.data(), projected_face_vertices.data());

    projected_face_vertices_2d.resize(2*n_face_vertices);
    for (int n = 0; n < n_face_vertices; ++n)
    {
        projected_face_vertices_2d[2*n+0] = projected_face_vertices[3*n+0];
        projected_face_vertices_2d[2*n+1] = projected_face_vertices[3*n+1];
    }

    /* DEBUG
    std::cout << std::endl;
    std::cout << "R: " << std::endl;
    for (int n = 0; n < 3; ++n)
    {
        std::cout << this->R[n+3*0] << ", " << this->R[n+3*1] << ", " << this->R[n+3*2] << std::endl;
    }

    std::cout << std::endl;
    for (int n = 0; n < n_face_vertices; ++n)
    {
        math::matrix::print_inline(3, &face_vertices[3*n]); std::cout << std::endl;
    }

    std::cout << std::endl;
    for (int n = 0; n < n_face_vertices; ++n)
    {
        math::matrix::print_inline(3, &projected_face_vertices[3*n]); std::cout << std::endl;
    }
    */
    // ----------------------------------------------------------------

    // GATHER THE BOUNDARY NODES --------------------------------------
    node_pick = 0;
    for (int fe = 0; fe < n_face_edges; ++fe)
    {
        ed = this->ed_conn[fe];
        const VoronoiEdge & edge = edges[ed];

        n_edge_elems = edge.m_conn.size();

        if (this->ed_ori[fe] == +1)
        {
            node_pick = 0;
        }
        else if (this->ed_ori[fe] == -1)
        {
            node_pick = 1;
        }
        else
        {
            io::error("mesh.cpp - VoronoiFace::build_mesh",
                      "Unexpected value of the edge orientation: "+std::to_string(this->ed_ori[fe]));
        }

        for (int n = 0; n < n_edge_elems; ++n)
        {
            e = edge.m_conn[n];
            elem_offset = global_mesh_offset[e];
            elem_conn = &global_mesh_conn[elem_offset];
            nodes.push_back(global_mesh_nodes[3*elem_conn[node_pick]+0]);
            nodes.push_back(global_mesh_nodes[3*elem_conn[node_pick]+1]);
            nodes.push_back(global_mesh_nodes[3*elem_conn[node_pick]+2]);
            bou_nodes.push_back(elem_conn[node_pick]);

        }
    }

    n_boundary_nodes = nodes.size()/3;

    /* DEBUG
    std::cout << "bou_nodes:" << std::endl;
    for (int n = 0; n < n_boundary_nodes; ++n)
    {
        std::cout << bou_nodes[n] << std::endl;
    }
    */

    // PROJECT THE BOUNDARY NODES ONTO THE FACE LOCAL REF SYSTEM
    projected_nodes.resize(nodes.size());
    math::matrix::multiply(3, 3, n_boundary_nodes, this->R, nodes.data(), projected_nodes.data());
    z = projected_nodes[3*0+2];

    projected_nodes_2d.resize(2*n_boundary_nodes);
    for (int n = 0; n < n_boundary_nodes; ++n)
    {
        projected_nodes_2d[2*n+0] = projected_nodes[3*n+0];
        projected_nodes_2d[2*n+1] = projected_nodes[3*n+1];
    }

    /* DEBUG
    std::cout << std::endl;
    for (int n = 0; n < n_boundary_nodes; ++n)
    {
        math::matrix::print_inline(3, &nodes[3*n]); std::cout << std::endl;
    }

    std::cout << std::endl;
    for (int n = 0; n < n_boundary_nodes; ++n)
    {
        math::matrix::print_inline(3, &projected_nodes[3*n]); std::cout << std::endl;
    }
    */
    // ----------------------------------------------------------------

    // MESH THE POLYGONAL AREA ----------------------------------------
    utils::mesh_polygon(mesh_size, projected_face_vertices_2d,
                        projected_nodes_2d,
                        nodes_2d, conn);

    n_nodes = nodes_2d.size()/2;
    projected_nodes.resize(3*n_nodes);

    for (int n = 0; n < n_nodes; ++n)
    {
        projected_nodes[3*n+0] = nodes_2d[2*n+0];
        projected_nodes[3*n+1] = nodes_2d[2*n+1];
        projected_nodes[3*n+2] = z;
    }
    // ----------------------------------------------------------------

    // PROJECT THE NODES BACK TO THE GLOBAL REFERENCE SYSTEM ----------
    nodes.resize(projected_nodes.size());
    math::matrix::transpose(3, 3, this->R, RT);
    math::matrix::multiply(3, 3, n_nodes, RT, projected_nodes.data(), nodes.data());
    // ----------------------------------------------------------------
}
// ====================================================================
// ####################################################################



// VORONOI CELL CLASS #################################################
// CONSTRUCTOR ========================================================
VoronoiCell::VoronoiCell()
:
state(0)
{

}
// ====================================================================

// BUILD MESH =========================================================
void VoronoiCell::build_mesh(const double mesh_size,
                             const std::vector<VoronoiVertex> & vertices,
                             const std::vector<VoronoiEdge> & edges,
                             const std::vector<VoronoiFace> & faces,
                             const std::vector<double> & global_mesh_nodes,
                             const std::vector<int> & global_mesh_offset,
                             const std::vector<int> & global_mesh_conn,
                             std::vector<double> & nodes, std::vector<int> & conn,
                             std::vector<int> & bou_nodes) const
{
    // PARAMETERS -----------------------------------------------------
    const int n_cell_faces = this->f_conn.size();
    // ----------------------------------------------------------------

    // VARIABLES ------------------------------------------------------
    std::vector<std::array<double, 4>> planes;
    int n_face_elems, e, elem_offset;
    const int * elem_conn;
    std::vector<int>::iterator iter;
    std::vector<double> boundary_nodes;
    std::vector<int> boundary_conn;
    int n_boundary_nodes, n_boundary_triangles;
    // ----------------------------------------------------------------

    // CLEAR MEMORY -------
    nodes.clear();
    conn.clear();
    bou_nodes.clear();
    boundary_nodes.clear();
    boundary_conn.clear();
    // --------------------

    // GATHER THE PARAMETERS OF THE PLANES DEFINING THE VORONOI CELL --
    for (int cf = 0; cf < n_cell_faces; ++cf)
    {
        const int f = this->f_conn[cf];
        const VoronoiFace & face = faces[f];
        const int v = face.v_conn[0];
        const double * V = vertices[v].x;

        std::array<double, 4> plane = {face.R[2+3*0], face.R[2+3*1], face.R[2+3*2], 0.0};
        plane[3] = -(plane[0]*V[0]+plane[1]*V[1]+plane[2]*V[2]);

        if (this->f_ori[cf] < 0)
        {
            plane[0] = -plane[0];
            plane[1] = -plane[1];
            plane[2] = -plane[2];
            plane[3] = -plane[3];
        }

        planes.push_back(plane);
    }
    // ----------------------------------------------------------------

    // GATHER THE BOUNDARY NODES --------------------------------------
    for (int cf = 0; cf < n_cell_faces; ++cf)
    {
        const int f = this->f_conn[cf];
        const VoronoiFace & face = faces[f];
        n_face_elems = face.m_conn.size();

        for (int n = 0; n < n_face_elems; ++n)
        {
            e = face.m_conn[n];
            elem_offset = global_mesh_offset[e];
            elem_conn = &global_mesh_conn[elem_offset];
            bou_nodes.push_back(elem_conn[0]);
            bou_nodes.push_back(elem_conn[1]);
            bou_nodes.push_back(elem_conn[2]);
            boundary_conn.push_back(elem_conn[0]);
            boundary_conn.push_back(elem_conn[1]);
            boundary_conn.push_back(elem_conn[2]);
        }
    }
    std::sort(bou_nodes.begin(), bou_nodes.end());
    iter = std::unique(bou_nodes.begin(), bou_nodes.end());
    bou_nodes.erase(iter, bou_nodes.end());

    n_boundary_nodes = bou_nodes.size();
    for (int n = 0; n < n_boundary_nodes; ++n)
    {
        boundary_nodes.push_back(global_mesh_nodes[3*bou_nodes[n]+0]);
        boundary_nodes.push_back(global_mesh_nodes[3*bou_nodes[n]+1]);
        boundary_nodes.push_back(global_mesh_nodes[3*bou_nodes[n]+2]);
    }

    n_boundary_triangles = boundary_conn.size()/3;
    for (int n = 0; n < 3*n_boundary_triangles; ++n)
    {
        iter = std::find(bou_nodes.begin(), bou_nodes.end(), boundary_conn[n]);
        boundary_conn[n] = std::distance(bou_nodes.begin(), iter);
    }
    // ----------------------------------------------------------------

    // MESH THE POLYHEDRAL DOMAIN -------------------------------------
    utils::c_mesh_convex_polyhedron(this->centroid,
                                    boundary_nodes, boundary_conn,
                                    nodes, conn);
    // ----------------------------------------------------------------

    /* DEBUG
    std::cout << std::endl;
    for (int n = 0; n < n_boundary_nodes; ++n)
    {
        math::matrix::print_inline(3, &boundary_nodes[3*n]); std::cout << std::endl;
    }
    */
    /* DEBUG
    std::cout << std::endl;
    for (int n = 0; n < n_boundary_triangles; ++n)
    {
        math::matrix::print_inline_int(3, &boundary_conn[3*n]); std::cout << std::endl;
    }
    */
}
// ====================================================================
// ####################################################################



// MESH CLASS #########################################################
// CONSTRUCTOR ========================================================
Mesh::Mesh()
:
state(0),
tol(0.0)
{

}
// ====================================================================

// DESTRUCTOR =========================================================
Mesh::~Mesh()
{
    this->free();
}
// ====================================================================

// FREE MEMORY ========================================================
void Mesh::clear_geometry()
{
    this->vertices.clear();
    this->edges.clear();
    this->faces.clear();
    this->cells.clear();
}

void Mesh::clear_mesh()
{
    this->nodes.clear();
    this->etype.clear();
    this->ghost.clear();
    this->offset.clear();
    this->conn.clear();
    this->parents.clear();
}

void Mesh::free()
{
    this->clear_mesh();

    this->clear_geometry();

    this->tol = 0.0;

    this->seeds.clear();

    const int n_vc = this->vc_ptrs.size();
    for (int c = 0; c < n_vc; ++c)
    if (this->vc_ptrs[c] != nullptr)
    {
        delete this->vc_ptrs[c];
    }
    this->vc_ptrs.clear();
}
// ====================================================================

// INITIALIZATION: TOLERANCE ==========================================
void Mesh::init_tolerance(const double rel_tol)
{
    // PARAMETERS ---------------------------
    const int n_cells = this->vc_ptrs.size();
    // --------------------------------------

    // VARIABLES ---------------------------------------------
    const double * A, * B;
    double shortest_edge = std::numeric_limits<double>::max();
    voro::voronoicell_neighbor * vc;
    std::vector<double> ver;
    // -------------------------------------------------------
            
    // LOOP OVER THE VORO++ NEIGHBOR CELL POINTERS --------------------
    for (int c = 0; c < n_cells; ++c)
    {
        if (this->vc_ptrs[c] != nullptr)
        {
            vc = this->vc_ptrs[c];
            vc->vertices(0.0, 0.0, 0.0, ver);

            for (int i = 0; i < vc->p; ++i)
            for (int j = 0; j < vc->nu[i]; ++j)
            {
                if ((vc->ed[i][j]) > i)
                {

                }
                else
                {
                    A = &ver[3*i];
                    B = &ver[3*(vc->ed[i][j])];
                    shortest_edge = std::min(shortest_edge, math::L2_distance(A, B, 3));
                }
            }
        }
    }
    // ----------------------------------------------------------------

    this->tol = shortest_edge*rel_tol;
}
// ====================================================================

// INITIALIZATION: VORONOI CELLS ======================================
void Mesh::init_voronoi_cells(const std::vector<int> & ids)
{
    // CLEAR PREVIOUS DATA
    this->cells.clear();
    // -------------------

    // PARAMETERS ---------------------------
    const int n_cells = this->vc_ptrs.size();
    // --------------------------------------

    // VARIABLES -------------------
    voro::voronoicell_neighbor * vc;
    // -----------------------------

    // LOOP OVER THE VORO++ NEIGHBOR CELL POINTERS --------------------
    for (int c = 0; c < n_cells; ++c)
    {
        vc = this->vc_ptrs[c];

        // CREATE A VORONOI CELL
        VoronoiCell cell;

        if (this->vc_ptrs[c] != nullptr)
        {
            // Set the id
            cell.id = ids[c];

            // Set the neighbors
            vc->neighbors(cell.nbr);

            // Compute the centroid
            vc->centroid(cell.centroid[0], cell.centroid[1], cell.centroid[2]);
            cell.centroid[0] += this->seeds[3*c+0];
            cell.centroid[1] += this->seeds[3*c+1];
            cell.centroid[2] += this->seeds[3*c+2];

        }
        else
        {
            // Set the id
            cell.id = -1;

            // Set the neighbors
            cell.nbr.clear();

            // Compute the centroid
            cell.centroid[0] = 0.0;
            cell.centroid[1] = 0.0;
            cell.centroid[2] = 0.0;
        }

        // Push the cell
        this->cells.push_back(cell);
    }
    // ----------------------------------------------------------------

    // MAKE SURE THE CONNECTIVITY IS CONSISTENT -----------------------
    for (int c = 0; c < n_cells; ++c)
    {
        VoronoiCell & cell = this->cells[c];
        const int n_nbr = cell.nbr.size();
            
        for (int n = 0; n < n_nbr; ++n)
        {
            int nbr_id = cell.nbr[n];

            // WALL
            if (nbr_id < 0)
            {

            }
            // THE CELL IS NEIGHBOR WITH ITSELF
            else if (nbr_id == cell.id)
            {
                io::error("mesh.cpp - Mesh::init_voronoi_cells",
                          "Cells neighboring with themselves must be handled yet.");
            }
            // AN ACTUAL NEIGHBORING CELL
            else
            {
                int nbr_c = 0;
                bool found = false;
                while ((!found) && (nbr_c < n_cells))
                {
                    const VoronoiCell & nbr_cell = this->cells[nbr_c];
                    if (nbr_cell.id == nbr_id)
                    {
                        found = true;
                        cell.nbr[n] = nbr_c;
                    }
                    else
                    {
                        nbr_c += 1;
                    }
                }

                if (!found)
                {
                    io::error("mesh.cpp - Mesh::init_voronoi_cells",
                              "The neighbor cells information is not consistent.");
                }
            }
        }
    }
    // ----------------------------------------------------------------
}
// ====================================================================

// INITIALIZATION: VORONOI FACES ======================================
void Mesh::init_voronoi_faces()
{
    // CLEAR PREVIOUS DATA
    this->faces.clear();
    // -------------------

    // PARAMETERS -------------------------
    const int n_cells = this->cells.size();
    // ------------------------------------

    // VARIABLES --------------------
    voro::voronoicell_neighbor * vc;
    std::vector<double> un, ver;
    std::vector<int> face_vertices;
    int face_vertices_pos;
    // ------------------------------

    // LOOP OVER THE VORONOI CELLS ------------------------------------
    for (int c = 0; c < n_cells; ++c)
    {
        if (this->vc_ptrs[c] != nullptr)
        {
            vc = this->vc_ptrs[c];
            vc->normals(un);
            vc->vertices(0.0, 0.0, 0.0, ver);
            vc->face_vertices(face_vertices);

            VoronoiCell & cell = this->cells[c];
            const int n_nbr = cell.nbr.size();
            const int n_cell_faces = n_nbr;

            cell.f_conn.resize(n_cell_faces);
            cell.f_ori.resize(n_cell_faces);

            face_vertices_pos = 0;
            for (int f = 0; f < n_cell_faces; ++f)
            {
                const int nbr_c = cell.nbr[f];
                const bool face_is_wall = (nbr_c < 0);
                const bool new_face = (face_is_wall || (nbr_c > c));

                // ADD A NEW FACE
                if (new_face)
                {
                    // CREATE A FACE
                    VoronoiFace face;

                    // Set the face type (as a boundary)
                    face.bou_type = ((face_is_wall)? BOU_TYPE_WALL : BOU_TYPE_INTERFACE);

                    // Set the parent cells
                    face.parent_cells[0] = c;
                    face.parent_cells[1] = nbr_c;

                    // Eval a local reference system
                    face.eval_RS(ver.data(), face_vertices[face_vertices_pos], &face_vertices[face_vertices_pos+1]);

                    // Store the index in the cell's face connectivity
                    cell.f_conn[f] = this->faces.size();

                    // Store the orientation of the face
                    const double * V0 = &ver[3*face_vertices[face_vertices_pos+1]];
                    const double * CC = cell.centroid;
                    const double CV[3] = {V0[0]-CC[0], V0[1]-CC[1], V0[2]-CC[2]};
                    const double fn[3] = {face.R[2+3*0], face.R[2+3*1], face.R[2+3*2]};
                    const double tmp = CV[0]*fn[0]+CV[1]*fn[1]+CV[2]*fn[2];
                    cell.f_ori[f] = (tmp > 0.0) ? +1 : -1;

                    // Push the face
                    this->faces.push_back(face);
                }
                // THE FACE ALREADY EXISTS
                else
                {
                    // FIND THE FACE
                    const int n_faces = this->faces.size();

                    bool found = false;
                    int ff = 0;
                    while ((!found) && (ff < n_faces))
                    {
                        const VoronoiFace & face = this->faces[ff];

                        if ((face.parent_cells[0] == nbr_c) &&
                            (face.parent_cells[1] == c) &&
                            (face.bou_type != BOU_TYPE_WALL))
                        {
                            found = true;

                            // Store the index in the cell's face
                            // connectivity
                            cell.f_conn[f] = ff;

                            // Store the orientation of the face
                        }
                        else
                        {
                            ff += 1;
                        }
                    }

                    if (!found)
                    {
                        io::error("mesh.cpp - Mesh::init_voronoi_faces",
                                "The neighbor cells information is not consistent.");
                    }
                }

                // MOVE TO THE NEXT FACE
                face_vertices_pos += face_vertices[face_vertices_pos]+1;
            }
        }
    }
    // ----------------------------------------------------------------
}
// ====================================================================

// BUILD MESH: BASE ===================================================
void Mesh::build(const int dim, const double mesh_size)
{
    this->build_bottom_up(dim, mesh_size);
}
// ====================================================================

// BUILD MESH: BOTTOM-UP ==============================================
void Mesh::build_bottom_up(const int dim, const double mesh_size)
{
    // PARAMETERS -----------------------------------------------------
    const int n_cells = this->cells.size();
    const int n_faces = this->faces.size();
    // ----------------------------------------------------------------

    // VARIABLES ------------------------------------------------------
    const double * S, * A, * B;
    voro::voronoicell_neighbor * vc;
    std::vector<double> ver;
    bool foundA, foundB, found;
    int vA, vB, ed, n_vertices, n_edges;

    int nodes_pos, aux_offset;
    // ----------------------------------------------------------------

    // CLEAR PREVIOUS DATA
    this->vertices.clear();
    this->edges.clear();
    // --------------------

    // CONSTRUCT THE UNIQUE LIST OF VERTICES AND EDGES ----------------
    for (int c = 0; c < n_cells; ++c)
    {
        if (this->vc_ptrs[c] != nullptr)
        {
            vc = this->vc_ptrs[c];
            S = &this->seeds[3*c];

            vc->vertices(S[0], S[1], S[2], ver);

            VoronoiCell & cell = this->cells[c];
            cell.v_conn.resize(ver.size()/3);

            for (int i = 0; i < vc->p; ++i)
            for (int j = 0; j < vc->nu[i]; ++j)
            {
                if ((vc->ed[i][j]) > i)
                {

                }
                else
                {
                    A = &ver[3*i];
                    B = &ver[3*(vc->ed[i][j])];

                    // CHECK WHETHER WE HAVE ALREADY ADDED THE VERTEX A
                    n_vertices = this->vertices.size();
                    foundA = false;
                    vA = 0;
                    while ((!foundA) && (vA < n_vertices))
                    {
                        const VoronoiVertex & vertex = this->vertices[vA];

                        if (math::L2_distance(A, vertex.x, 3) < this->tol)
                        {
                            foundA = true;
                        }
                        else
                        {
                            vA += 1;
                        }
                    }

                    if (!foundA)
                    {
                        VoronoiVertex vertex;
                        vertex.x[0] = A[0];
                        vertex.x[1] = A[1];
                        vertex.x[2] = A[2];
                        this->vertices.push_back(vertex);
                    }

                    cell.v_conn[i] = vA;

                    // CHECK WHETHER WE HAVE ALREADY ADDED THE VERTEX B
                    n_vertices = this->vertices.size();
                    foundB = false;
                    vB = 0;
                    while ((!foundB) && (vB < n_vertices))
                    {
                        const VoronoiVertex & vertex = this->vertices[vB];

                        if (math::L2_distance(B, vertex.x, 3) < this->tol)
                        {
                            foundB = true;
                        }
                        else
                        {
                            vB += 1;
                        }
                    }

                    if (!foundB)
                    {
                        VoronoiVertex vertex;
                        vertex.x[0] = B[0];
                        vertex.x[1] = B[1];
                        vertex.x[2] = B[2];
                        this->vertices.push_back(vertex);
                    }

                    cell.v_conn[vc->ed[i][j]] = vB;

                    // CHECK WHETHER WE HAVE ALREADY ADDED THE EDGE AB
                    n_edges = this->edges.size();
                    found = false;
                    ed = 0;
                    while ((!found) && (ed < n_edges))
                    {
                        const VoronoiEdge & edge = this->edges[ed];

                        if ((edge.v_conn[0] == vA) && (edge.v_conn[1] == vB))
                        {
                            found = true;
                            cell.ed_conn.push_back(ed);
                            cell.ed_ori.push_back(+1);
                        }
                        else if ((edge.v_conn[0] == vB) && (edge.v_conn[1] == vA))
                        {
                            found = true;
                            cell.ed_conn.push_back(ed);
                            cell.ed_ori.push_back(-1);
                        }
                        else
                        {
                            ed += 1;
                        }
                    }

                    if (!found)
                    {
                        VoronoiEdge edge;
                        edge.v_conn[0] = vA;
                        edge.v_conn[1] = vB;
                        cell.ed_conn.push_back(this->edges.size());
                        cell.ed_ori.push_back(+1);
                        this->edges.push_back(edge);
                    }
                }
            }
        }
    }

    // NUMBER OF VERTICES AND EDGES
    n_vertices = this->vertices.size();
    n_edges = this->edges.size();
    // ----------------------------------------------------------------

    // CONSTRUCT THE VERTEX AND EDGE CONNECTIVITY OF THE FACES --------
    for (int c = 0; c < n_cells; ++c)
    {
        if (this->vc_ptrs[c] != nullptr)
        {
            vc = this->vc_ptrs[c];
            std::vector<int> face_vertices;
            int face_vertices_pos;
            vc->face_vertices(face_vertices);

            const VoronoiCell & cell = this->cells[c];
            const int n_cell_faces = cell.f_conn.size();

            face_vertices_pos = 0;
            for (int cf = 0; cf < n_cell_faces; ++cf)
            {
                if ((cell.nbr[cf] > c) || (cell.nbr[cf] < 0))
                {
                    VoronoiFace & face = this->faces[cell.f_conn[cf]];

                    // SET THE VERTEX CONNECTIVITY
                    const int n_face_vertices = face_vertices[face_vertices_pos];
                    face.v_conn.resize(n_face_vertices);
                    for (int v = 0; v < n_face_vertices; ++v)
                    {
                        const int vv = face_vertices[face_vertices_pos+1+v];
                        face.v_conn[v] = cell.v_conn[vv];
                    }
                    
                    // SET THE EDGE CONNECTIVITY
                    for (int v = 0; v < n_face_vertices; ++v)
                    {
                        vA = face.v_conn[v];
                        vB = (v < (n_face_vertices-1))? face.v_conn[v+1] : face.v_conn[0];

                        found = false;
                        ed = 0;
                        while ((!found) && (ed < n_edges))
                        {
                            const VoronoiEdge & edge = this->edges[ed];

                            if ((edge.v_conn[0] == vA) && (edge.v_conn[1] == vB))
                            {
                                found = true;
                                face.ed_conn.push_back(ed);
                                face.ed_ori.push_back(+1);
                            }
                            else if ((edge.v_conn[0] == vB) && (edge.v_conn[1] == vA))
                            {
                                found = true;
                                face.ed_conn.push_back(ed);
                                face.ed_ori.push_back(-1);
                            }
                            else
                            {
                                ed += 1;
                            }
                        }

                        if (!found)
                        {
                            io::error("mesh.cpp - Mesh::export_vtk",
                                    "Could not find a face' edge.");
                        }
                    }
                }

                // MOVE TO THE NEXT FACE
                face_vertices_pos += face_vertices[face_vertices_pos]+1;
            }
        }
    }
    // ----------------------------------------------------------------

    // CLEAR MESH MEMORY
    this->clear_mesh();
    // ----------------

    // ADD THE VERTICES' NODES AND ELEMENTS ---------------------------
    // MAKE ROOM
    this->nodes.resize(3*n_vertices);
    this->etype.resize(n_vertices);
    this->ghost.resize(n_vertices);
    this->offset.resize(n_vertices);
    this->conn.resize(n_vertices);
    this->parents.resize(n_vertices);
    
    // ADD NODES AND ELEMENTS
    for (int v = 0; v < n_vertices; ++v)
    {
        VoronoiVertex & vertex = this->vertices[v];

        // Nodes
        this->nodes[3*v+0] = vertex.x[0];
        this->nodes[3*v+1] = vertex.x[1];
        this->nodes[3*v+2] = vertex.x[2];

        // Element type
        this->etype[v] = VTK_VERTEX;

        // Ghost type
        this->ghost[v] = 0;

        // Offset
        this->offset[v] = v;

        // Connectivity
        this->conn[v] = v;

        // Parent geometrical entity
        this->parents[v] = v;

        // Update vertex info
        vertex.m_conn = v;
    }

    if (dim < 1) return;
    // ----------------------------------------------------------------

    // ADD THE EDGES' NODES AND ELEMENTS ------------------------------
    // AUXILIARY NODES OFFSET POSITION
    nodes_pos = n_vertices;
    aux_offset = n_vertices;
    
    for (int ed = 0; ed < n_edges; ++ed)
    {
        VoronoiEdge & edge = this->edges[ed];
        std::vector<double> edge_nodes;
        std::vector<int> edge_conn;
        const int n_elems = this->etype.size();

        // Build the mesh of the edge
        edge.build_mesh(mesh_size,
                        this->vertices,
                        edge_nodes, edge_conn);

        // Some data
        const int n_edge_nodes = edge_nodes.size()/3;
        const int n_edge_elems = edge_conn.size()/2;

        // Nodes
        for (int n = 2; n < n_edge_nodes; ++n)
        {
            this->nodes.push_back(edge_nodes[3*n+0]);
            this->nodes.push_back(edge_nodes[3*n+1]);
            this->nodes.push_back(edge_nodes[3*n+2]);
        }

        for (int e = 0; e < n_edge_elems; ++e)
        {
            // Element type
            this->etype.push_back(VTK_LINE);

            // Ghost type
            this->ghost.push_back(0);

            // Offset
            this->offset.push_back(aux_offset);
            aux_offset += 2;

            // Parent geometrical entity
            this->parents.push_back(ed);
        }

        // Modify the connectivity of the edge
        const int n_old_nodes = 2;
        for (int e = 0; e < n_edge_elems; ++e)
        {
            edge_conn[2*e+0] += nodes_pos-n_old_nodes;
            edge_conn[2*e+1] += nodes_pos-n_old_nodes;
        }
        edge_conn[2*0+0] = this->vertices[edge.v_conn[0]].m_conn;
        edge_conn[2*(n_edge_elems-1)+1] = this->vertices[edge.v_conn[1]].m_conn;

        nodes_pos += n_edge_elems-1;

        // Connectivity
        for (int e = 0; e < n_edge_elems; ++e)
        {
            this->conn.push_back(edge_conn[2*e+0]);
            this->conn.push_back(edge_conn[2*e+1]);
        }

        // Update edge info
        edge.m_conn.resize(n_edge_elems);
        for (int e = 0; e < n_edge_elems; ++e)
        {
            edge.m_conn[e] = n_elems+e;
        }
    }

    if (dim < 2) return;
    // ----------------------------------------------------------------

    // ADD THE FACES' NODES AND ELEMENTS ------------------------------
    for (int f = 0; f < n_faces; ++f)
    {
        VoronoiFace & face = this->faces[f];
        std::vector<double> face_nodes;
        std::vector<int> face_conn, bou_nodes;
        const int n_elems = this->etype.size();
        int n_face_nodes, face_conn_size, n_face_elems, n_boundary_nodes;
        int n_new_nodes;

        // Build the mesh of the face
        face.build_mesh(mesh_size,
                        this->vertices, this->edges,
                        this->nodes, this->offset, this->conn,
                        face_nodes, face_conn, bou_nodes);

        n_face_nodes = face_nodes.size()/3;
        face_conn_size = face_conn.size();
        n_face_elems = face_conn_size/3;
        n_boundary_nodes = bou_nodes.size();
        n_new_nodes = n_face_nodes-n_boundary_nodes;

        // Nodes
        for (int n = n_boundary_nodes; n < n_face_nodes; ++n)
        {
            this->nodes.push_back(face_nodes[3*n+0]);
            this->nodes.push_back(face_nodes[3*n+1]);
            this->nodes.push_back(face_nodes[3*n+2]);
        }

        for (int e = 0; e < n_face_elems; ++e)
        {
            // Element type
            this->etype.push_back(VTK_TRIANGLE);

            // Ghost type
            this->ghost.push_back(0);

            // Offset
            this->offset.push_back(aux_offset);
            aux_offset += 3;

            // Parent geometrical entity
            this->parents.push_back(f);
        }

        /* DEBUG
        std::cout << "face_conn:" << std::endl;
        for (int n = 0; n < n_face_elems; ++n)
        {
            math::matrix::print_inline_int(3, &face_conn[3*n]); std::cout << std::endl;
        }
        */

        // Modify the connectivity of the face
        for (int c = 0; c < face_conn_size; ++c)
        {
            if (face_conn[c] < n_boundary_nodes)
            {
                face_conn[c] = bou_nodes[face_conn[c]];
            }
            else
            {
                face_conn[c] = face_conn[c]-n_boundary_nodes+nodes_pos;
            }
        }

        nodes_pos += n_new_nodes;

        /* DEBUG
        std::cout << "face_nodes:" << std::endl;
        for (int n = 0; n < n_face_nodes; ++n)
        {
            math::matrix::print_inline(3, &face_nodes[3*n]); std::cout << std::endl;
        }
        std::cout << "face_conn:" << std::endl;
        for (int n = 0; n < n_face_triangles; ++n)
        {
            math::matrix::print_inline_int(3, &face_conn[3*n]); std::cout << std::endl;
        }
        std::cout << "bou_nodes:" << std::endl;
        for (int n = 0; n < n_boundary_nodes; ++n)
        {
            std::cout << bou_nodes[n] << std::endl;
        }
        std::cout << std::endl;
        */

        // Connectivity
        for (int e = 0; e < n_face_elems; ++e)
        {
            this->conn.push_back(face_conn[3*e+0]);
            this->conn.push_back(face_conn[3*e+1]);
            this->conn.push_back(face_conn[3*e+2]);
        }

        // Update face info
        face.m_conn.resize(n_face_elems);
        for (int e = 0; e < n_face_elems; ++e)
        {
            face.m_conn[e] = n_elems+e;
        }
    }

    if (dim < 3) return;
    // ----------------------------------------------------------------


    // ADD THE POLYHEDRA' NODES AND ELEMENTS --------------------------
    for (int c = 0; c < n_cells; ++c)
    {
        VoronoiCell & cell = this->cells[c];
        std::vector<double> cell_nodes;
        std::vector<int> cell_conn, bou_nodes;
        const int n_elems = this->etype.size();
        int n_cell_nodes, cell_conn_size, n_cell_elems, n_boundary_nodes;
        int n_new_nodes;

        // Build the mesh of the cell
        cell.build_mesh(mesh_size,
                        this->vertices, this->edges, this->faces,
                        this->nodes, this->offset, this->conn,
                        cell_nodes, cell_conn, bou_nodes);

        n_cell_nodes = cell_nodes.size()/3;
        cell_conn_size = cell_conn.size();
        n_cell_elems = cell_conn_size/4;
        n_boundary_nodes = bou_nodes.size();
        n_new_nodes = n_cell_nodes-n_boundary_nodes;

        // Nodes
        for (int n = n_boundary_nodes; n < n_cell_nodes; ++n)
        {
            this->nodes.push_back(cell_nodes[3*n+0]);
            this->nodes.push_back(cell_nodes[3*n+1]);
            this->nodes.push_back(cell_nodes[3*n+2]);
        }

        for (int e = 0; e < n_cell_elems; ++e)
        {
            // Element type
            this->etype.push_back(VTK_TETRA);

            // Ghost type
            this->ghost.push_back(0);

            // Offset
            this->offset.push_back(aux_offset);
            aux_offset += 4;

            // Parent geometrical entity
            this->parents.push_back(c);
        }

        // Modify the connectivity of the cell
        for (int c = 0; c < cell_conn_size; ++c)
        {
            if (cell_conn[c] < n_boundary_nodes)
            {
                cell_conn[c] = bou_nodes[cell_conn[c]];
            }
            else
            {
                cell_conn[c] = cell_conn[c]-n_boundary_nodes+nodes_pos;
            }
        }

        nodes_pos += n_new_nodes;

        // Connectivity
        for (int e = 0; e < n_cell_elems; ++e)
        {
            this->conn.push_back(cell_conn[4*e+0]);
            this->conn.push_back(cell_conn[4*e+1]);
            this->conn.push_back(cell_conn[4*e+2]);
            this->conn.push_back(cell_conn[4*e+3]);
        }

        // Update face info
        cell.m_conn.resize(n_cell_elems);
        for (int e = 0; e < n_cell_elems; ++e)
        {
            cell.m_conn[e] = n_elems+e;
        }
    }
    // ----------------------------------------------------------------
}
// ====================================================================

// BUILD MESH: CELL-BASED =============================================
void Mesh::build_cell_based(const int dim, const double mesh_size)
{

}
// ====================================================================

// EXPORT TO VTK FORMAT ===============================================
void Mesh::export_vtk(const std::string & filepath)
{
    // PARAMETERS -----------------------------------------------------
    const vtk::Cell_conn_t n_VTK_nodes = this->nodes.size()/3;
    const vtk::Cell_offs_t n_VTK_cells = this->etype.size();

    // CELL DATA
    std::vector<std::string> VTK_cell_fields_name = 
    {
        "DIM",
        "PARENT"
    };
    const int n_cell_fields = VTK_cell_fields_name.size();
    // ----------------------------------------------------------------

    // VARIABLES ------------------------------------------------------
    long VTK_cell_conn_len;

    vtk::Cell_offs_t VTK_cells_pos;

    std::vector<vtk::Float_t> VTK_nodes;
    std::vector<vtk::Cell_conn_t> VTK_cell_conn;
    std::vector<vtk::Cell_offs_t> VTK_cell_offset;
    std::vector<vtk::Cell_type_t> VTK_cell_type;
    std::vector<std::vector<vtk::Int_t>> VTK_cell_fields(n_cell_fields);
    // ----------------------------------------------------------------

    // EVAL THE SIZE OF THE VTK CELL CONNECTIVITY LIST ----------------
    VTK_cell_conn_len = 0L;

    for (vtk::Cell_offs_t e = 0; e < n_VTK_cells; ++e)
    {
        if (this->etype[e] == VTK_VERTEX)
        {
            VTK_cell_conn_len += 1;
        }
        else if (this->etype[e] == VTK_LINE)
        {
            VTK_cell_conn_len += 2;
        }
        else if (this->etype[e] == VTK_TRIANGLE)
        {
            VTK_cell_conn_len += 3;
        }
        else if (this->etype[e] == VTK_TETRA)
        {
            VTK_cell_conn_len += 4;
        }
        else
        {
            io::error("mesh.cpp - Mesh::export_vtk",
                      "Unexpected element type: "+std::to_string(this->etype[e]));
        }
    }
    // ----------------------------------------------------------------

    // RESIZE MEMORY --------------------------------------------------
    VTK_nodes.resize(3*n_VTK_nodes);
    VTK_cell_conn.resize(VTK_cell_conn_len);
    VTK_cell_offset.resize(n_VTK_cells+1);
    VTK_cell_type.resize(n_VTK_cells);
    for (int f = 0; f < n_cell_fields; ++f)
    {
        VTK_cell_fields[f].resize(n_VTK_cells);
    }
    // ----------------------------------------------------------------

    // EVAL VTK NODES -------------------------------------------------
    for (vtk::Cell_conn_t n = 0; n < n_VTK_nodes; ++n)
    {
        VTK_nodes[3*n+0] = this->nodes[3*n+0];
        VTK_nodes[3*n+1] = this->nodes[3*n+1];
        VTK_nodes[3*n+2] = this->nodes[3*n+2];
    }
    // ----------------------------------------------------------------

    // EVAL VTK CELLS -------------------------------------------------
    VTK_cell_offset[0] = 0;
    
    VTK_cells_pos = 0;
    for (vtk::Cell_offs_t e = 0; e < n_VTK_cells; ++e)
    {
        if (this->etype[e] == VTK_VERTEX)
        {
            VTK_cell_conn[VTK_cells_pos] = this->conn[VTK_cells_pos];
            VTK_cells_pos += 1;
        }
        else if (this->etype[e] == VTK_LINE)
        {
            VTK_cell_conn[VTK_cells_pos+0] = this->conn[VTK_cells_pos+0];
            VTK_cell_conn[VTK_cells_pos+1] = this->conn[VTK_cells_pos+1];
            VTK_cells_pos += 2;
        }
        else if (this->etype[e] == VTK_TRIANGLE)
        {
            VTK_cell_conn[VTK_cells_pos+0] = this->conn[VTK_cells_pos+0];
            VTK_cell_conn[VTK_cells_pos+1] = this->conn[VTK_cells_pos+1];
            VTK_cell_conn[VTK_cells_pos+2] = this->conn[VTK_cells_pos+2];
            VTK_cells_pos += 3;
        }
        else if (this->etype[e] == VTK_TETRA)
        {
            VTK_cell_conn[VTK_cells_pos+0] = this->conn[VTK_cells_pos+0];
            VTK_cell_conn[VTK_cells_pos+1] = this->conn[VTK_cells_pos+1];
            VTK_cell_conn[VTK_cells_pos+2] = this->conn[VTK_cells_pos+2];
            VTK_cell_conn[VTK_cells_pos+3] = this->conn[VTK_cells_pos+3];
            VTK_cells_pos += 4;
        }
        else
        {
            io::error("mesh.cpp - Mesh::export_vtk",
                      "Unexpected element type: "+std::to_string(this->etype[e]));
        }

        VTK_cell_offset[e+1] = VTK_cells_pos;
        VTK_cell_type[e] = this->etype[e];
    }
    // ----------------------------------------------------------------

    // VTK CELL FIELDS ------------------------------------------------
    for (vtk::Cell_offs_t e = 0; e < n_VTK_cells; ++e)
    {
        // DIM
        if (this->etype[e] == VTK_VERTEX)
        {
            VTK_cell_fields[0][e] = 0;
        }
        else if (this->etype[e] == VTK_LINE)
        {
            VTK_cell_fields[0][e] = 1;
        }
        else if (this->etype[e] == VTK_TRIANGLE)
        {
            VTK_cell_fields[0][e] = 2;
        }
        else if (this->etype[e] == VTK_TETRA)
        {
            VTK_cell_fields[0][e] = 3;
        }
        else
        {
            io::error("mesh.cpp - Mesh::export_vtk",
                      "Unexpected element type: "+std::to_string(this->etype[e]));
        }

        // PARENT
        VTK_cell_fields[1][e] = this->parents[e];
    }
    // ----------------------------------------------------------------

    // DUMP DATA ------------------------------------------------------
    vtk::print_unstructured_data(filepath,
                                 n_VTK_nodes,
                                 n_VTK_cells,
                                 VTK_nodes,
                                 VTK_cell_conn,
                                 VTK_cell_offset,
                                 VTK_cell_type,
                                 VTK_cell_fields,
                                 VTK_cell_fields_name,
                                 {},
                                 {},
                                 "ascii");
    // ----------------------------------------------------------------
}
// ====================================================================
// ####################################################################

}