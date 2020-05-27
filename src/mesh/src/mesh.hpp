// mesh.hpp

#ifndef VOROMESH_MESH_HPP_
#define VOROMESH_MESH_HPP_

#include <vector>
#include <voro++.hh>

#include "utils_io.hpp"
#include "utils_math.hpp"
#include "utils_voropp.hpp"

#define BOU_TYPE_UNDEFINED -1
#define BOU_TYPE_WALL 0
#define BOU_TYPE_INTERFACE 1
#define BOU_TYPE_PERIODIC 2

namespace voromesh
{
// VORONOI VERTEX CLASS ###############################################
struct VoronoiVertex
{
    // DATA MEMBERS ===================================================
    int state;

    double x[3];

    int m_conn;
    // ================================================================

    // CONSTRUCTOR ====================================================
    VoronoiVertex();
    // ================================================================
};
// ####################################################################



// VORONOI EDGE CLASS #################################################
struct VoronoiEdge
{
    // DATA MEMBERS ===================================================
    int state;

    int v_conn[2];

    std::vector<int> m_conn;
    // ================================================================

    // CONSTRUCTOR ====================================================
    VoronoiEdge();
    // ================================================================

    // BUILD MESH =====================================================
    void build_mesh(const double mesh_size,
                    const std::vector<VoronoiVertex> & vertices,
                    std::vector<double> & nodes, std::vector<int> & conn) const;
    // ================================================================
};
// ####################################################################



// VORONOI FACE CLASS #################################################
struct VoronoiFace
{
    // DATA MEMBERS ===================================================
    int state;

    int bou_type;
    int parent_cells[2];
    double R[9];

    std::vector<int> v_conn;
    std::vector<int> ed_conn, ed_ori;

    std::vector<int> m_conn;
    // ================================================================

    // CONSTRUCTOR ====================================================
    VoronoiFace();
    // ================================================================

    // EVAL LOCAL REFERENCE SYSTEM ====================================
    void eval_RS(const double * vertices, const int n_vertices, const int * conn);
    // ================================================================

    // BUILD MESH =====================================================
    void build_mesh(const double mesh_size,
                    const std::vector<VoronoiVertex> & vertices,
                    const std::vector<VoronoiEdge> & edges,
                    const std::vector<double> & global_mesh_nodes,
                    const std::vector<int> & global_mesh_offset,
                    const std::vector<int> & global_mesh_conn,
                    std::vector<double> & nodes, std::vector<int> & conn,
                    std::vector<int> & bou_nodes) const;
    // ================================================================
};
// ####################################################################



// VORONOI CELL CLASS #################################################
struct VoronoiCell
{
    // DATA MEMBERS ===================================================
    int state;

    int id;
    std::vector<int> nbr;
    double centroid[3];

    std::vector<int> f_conn, f_ori;
    std::vector<int> ed_conn, ed_ori;
    std::vector<int> v_conn;

    std::vector<int> m_conn;
    // ================================================================

    // CONSTRUCTOR ====================================================
    VoronoiCell();
    // ================================================================

    // BUILD MESH =====================================================
    void build_mesh(const double mesh_size,
                    const std::vector<VoronoiVertex> & vertices,
                    const std::vector<VoronoiEdge> & edges,
                    const std::vector<VoronoiFace> & faces,
                    const std::vector<double> & global_mesh_nodes,
                    const std::vector<int> & global_mesh_offset,
                    const std::vector<int> & global_mesh_conn,
                    std::vector<double> & nodes, std::vector<int> & conn,
                    std::vector<int> & bou_nodes) const;
    // ================================================================
};
// ####################################################################



// MESH CLASS #########################################################
struct Mesh
{
    // DATA MEMBERS ===================================================
    int state;
    std::vector<voro::voronoicell_neighbor *> vc_ptrs;

    std::vector<double> seeds;

    double tol;
    std::vector<VoronoiCell> cells;
    std::vector<VoronoiFace> faces;
    std::vector<VoronoiEdge> edges;
    std::vector<VoronoiVertex> vertices;

    std::vector<double> nodes;
    std::vector<int> etype;
    std::vector<int> ghost;
    std::vector<int> offset;
    std::vector<int> conn;
    std::vector<int> parents;
    // ================================================================

    // CONSTRUCTOR ====================================================
    Mesh();
    // ================================================================

    // DESTRUCTOR =====================================================
    ~Mesh();
    // ================================================================

    // FREE MEMORY ====================================================
    void clear_geometry();
    void clear_mesh();
    void free();
    // ================================================================

    // INITIALIZATION FROM A VORO++ CONTAINER =========================
    template <typename CONTAINER>
    void init(CONTAINER & con)
    {
        // PARAMETERS -----------------------
        const int np = con.total_particles();
        // ----------------------------------

        // FREE MEMORY
        this->free();
        // -----------
        
        // INIT MEMORY ----------
        this->vc_ptrs.resize(np);
        this->seeds.resize(3*np);
        // ----------------------

        // COMPUTE THE VORONOI CELLS ----------------------------------
        int k = 0;
        double * pp;
        std::vector<int> ids(np, -1);
        voro::voronoicell_neighbor * vc;
        voro::c_loop_all cl(con);

        if (cl.start()) do
        {
            vc = new voro::voronoicell_neighbor();
            
            if (con.compute_cell(*vc, cl))
            {
                this->vc_ptrs[k] = vc;

                pp = con.p[cl.ijk]+3*cl.q;
                this->seeds[3*k+0] = pp[0];
                this->seeds[3*k+1] = pp[1];
                this->seeds[3*k+2] = pp[2];

                // CHECK THE ID OF THE VORONOI CELL IS NOT NEGATIVE
                if (cl.pid() < 0)
                {
                    io::error("mesh.hpp - Mesh::init",
                              "Voronoi cells cannot have negative ids.");
                }
                ids[k] = cl.pid();
            }
            else
            {
                delete vc;
                this->vc_ptrs[k] = nullptr;
                ids[k] = -1;
            }

            k += 1;
        }
        while (cl.inc());
        // ------------------------------------------------------------

        // SET THE TOLERANCE --
        this->init_tolerance();
        // --------------------

        // ADD THE VORONOI CELLS -----
        this->init_voronoi_cells(ids);
        // ---------------------------

        // ADD THE VORONOI FACES --
        this->init_voronoi_faces();
        // ------------------------
    }
    // ================================================================

    // INITIALIZATION: TOLERANCE ======================================
    void init_tolerance(const double rel_tol = 1.0e-5);
    // ================================================================

    // INITIALIZATION: VORONOI CELLS ==================================
    void init_voronoi_cells(const std::vector<int> & ids);
    // ================================================================

    // INITIALIZATION: VORONOI FACES ==================================
    void init_voronoi_faces();
    // ================================================================

    // AUXILIARY METHODS ==============================================
    int get_new_id() const;
    int get_cell_index(const int id) const;
    // ================================================================

    // BUILD MESH: BASE ===============================================
    void build(const int dim, const double mesh_size);
    // ================================================================

    // BUILD MESH: BOTTOM-UP ==========================================
    void build_bottom_up(const int dim, const double mesh_size);
    // ================================================================

    // BUILD MESH: CELL-BASED =========================================
    void build_cell_based(const int dim, const double mesh_size);
    // ================================================================

    // BUILD MESH: FOR BOUNDARY ELEMENT METHODS =======================
    void build_for_BEM(const int dim, const double mesh_size);
    // ================================================================

    // CUT CELL =======================================================
    void cut_cell_by_point_and_plane(const int id, const double * P, const double * un);
    void cut_cell_by_vector(const int id, const double * V);
    void cut_cell_by_plane(const int id, const double * plane);
    // ================================================================

    // UPDATE MESH ====================================================
    // ================================================================

    // EXPORT TO VTK FORMAT ===========================================
    void export_vtk(const std::string & filepath);
    // ================================================================
};
// ####################################################################

}

#endif