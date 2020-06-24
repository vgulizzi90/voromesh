// mesh.hpp

#ifndef VOROMESH_MESH_HPP_
#define VOROMESH_MESH_HPP_

#include <vector>
#include <algorithm>

#include <voro++.hh>

#include "utils_io.hpp"
#include "utils_math.hpp"
#include "utils_voropp.hpp"

#define BOU_TYPE_UNDEFINED -1
#define BOU_TYPE_WALL 0
#define BOU_TYPE_INTERFACE 1
#define BOU_TYPE_PERIODIC 2
#define BOU_TYPE_CRACK 3

#define BOU_TYPE_BAD_NEIGHBOR -1010
#define WALL_ID_BAD_NEIGHBOR -1010

namespace voromesh
{
// VORONOI VERTEX CLASS ###############################################
struct VoronoiVertex
{
    // DATA MEMBERS ===================================================
    int state;

    double x[3];
    std::vector<int> parent_edges;
    std::vector<int> parent_faces;
    std::vector<int> parent_cells;

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
    std::vector<int> parent_faces;
    std::vector<int> parent_cells;

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
    double R[9];

    std::vector<int> v_conn;
    std::vector<int> ed_conn, ed_ori;
    int parent_cells[2];

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



// VORONOI WALL CLASS #################################################
struct VoronoiWall
{
    // DATA MEMBERS ===================================================
    int state;

    int id;

    std::vector<int> f_conn;
    std::vector<int> c_conn;
    // ================================================================

    // CONSTRUCTOR ====================================================
    VoronoiWall();
    // ================================================================
};
// ####################################################################



// VORONOI CRACK CLASS ################################################
struct VoronoiCrack
{
    // DATA MEMBERS ===================================================
    int state;

    int id;

    std::vector<int> f_conn;
    std::vector<std::array<int, 2>> c_conn;
    // ================================================================

    // CONSTRUCTOR ====================================================
    VoronoiCrack();
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
    std::vector<VoronoiWall> walls;
    std::vector<VoronoiCrack> cracks;

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

        // ADD THE VERTICES AND THE EDGES ------
        this->init_voronoi_vertices_and_edges();
        // -------------------------------------

        // ADD THE VORONOI WALLS --
        this->init_voronoi_walls();
        // ------------------------
    }
    // ================================================================

    // ADD VORO++ CELL ================================================
    template<class VC>
    void add_vc(const VC & vc, const int new_id)
    {
        // PARAMETERS -------------------------
        const int n_cells = this->cells.size();
        // ------------------------------------

        // VARIABLES --------------------------------------------------
        int c;
        voro::voronoicell_neighbor * new_vc;
        bool new_vc_is_valid;
        // ------------------------------------------------------------

        // CHECK THE ID OF THE CELL -----------------------------------
        if (new_id < 0)
        {
            io::error("mesh.hpp - Mesh::init",
                      "Voronoi cells cannot have negative ids.");
        }
        c = this->get_cell_index(new_id, 0);
        if (c != n_cells)
        {
            io::warning("mesh.hpp - Mesh::init",
                        "A voronoi cell with id = "+std::to_string(new_id)+" already exists.");
        }
        // ------------------------------------------------------------

        // CREATE A NEW VORO++ CELL -----------------------------------
        new_vc = new voro::voronoicell_neighbor();
        *new_vc = vc;
        new_vc_is_valid = (new_vc != nullptr);
        // ------------------------------------------------------------

        // UPDATE THE vc_ptr VECTOR --------
        if (new_vc_is_valid)
        {
            this->vc_ptrs.push_back(new_vc);
        }
        // ---------------------------------

        // UPDATE THE seeds VECTOR ----
        if (new_vc_is_valid)
        {
            this->seeds.push_back(0.0);
            this->seeds.push_back(0.0);
            this->seeds.push_back(0.0);
        }
        // ----------------------------

        // UPDATE TOLERANCE -------
        if (new_vc_is_valid)
        {
            this->init_tolerance();
        }
        // ------------------------

        // RECOMPUTE THE VORONOI CELLS AND FACES ---
        if (new_vc_is_valid)
        {
            std::vector<int> ids(n_cells+1);
            for (int i = 0; i < n_cells; ++i)
            {
                ids[i] = this->cells[i].id;
            }
            ids[n_cells] = new_id;
            this->init_voronoi_cells(ids);
            this->init_voronoi_faces();
            this->init_voronoi_vertices_and_edges();
            this->init_voronoi_walls();
        }
        // -----------------------------------------
    }
    // ================================================================

    // AUXILIARY METHODS ==============================================
    int get_new_id() const;
    int get_new_wall_id() const;
    int get_cell_index(const int id, const int verbosity = 0) const;
    int get_wall_index(const int wall_id, const int verbosity = 0) const;
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

    // INITIALIZATION: VORONOI VERTICES AND EDGES =====================
    void init_voronoi_vertices_and_edges();
    // ================================================================

    // INITIALIZATION: WALLS ==========================================
    void init_voronoi_walls();
    // ================================================================

    // BUILD MESH: BASE ===============================================
    void build(const int dim, const double mesh_size);
    // ================================================================

    // BUILD MESH: BOTTOM-UP ==========================================
    void build_bottom_up(const int dim, const double mesh_size);
    // ================================================================

    // CUT CELL =======================================================
    void cut_cell_by_vector(const int id, const double * V);
    void cut_cell_by_plane(const int id, const double * plane);
    // ================================================================

    // CUT BASE METHODS ===============================================
    template<class DISTANCE_FUNCTION>
    void split_edges(const DISTANCE_FUNCTION & phi)
    {
        // PARAMETERS -------------------------------------------------
        const int n_edges_0 = this->edges.size();
        const int n_it = 100;
        // ------------------------------------------------------------

        // VARIABLES --------------------------------------------------
        // ------------------------------------------------------------

        // LOOP OVER THE EDGES ----------------------------------------
        for (int ed = 0; ed < n_edges_0; ++ed)
        {
            VoronoiEdge & edge = this->edges[ed];
            const double * A = this->vertices[edge.v_conn[0]].x;
            const double * B = this->vertices[edge.v_conn[1]].x;
            const double phiA = phi.eval(A);
            const double phiB = phi.eval(B);
            
            if ((std::abs(phiA) < this->tol) || (std::abs(phiA) < this->tol))
            {
                continue;
            }
            if (phiA*phiB > 0.0)
            {
                continue;
            }

            // FIND THE INTERSECTION BETWEEN THE EDGE AND THE ZERO
            // LEVEL SET OF THE DISTANCE FUNCTION
            double t;
            double X[3];
            double phiX, grad_phi[3], dphidt;
            int it;

            t = phiA/(phiA-phiB);
            it = 0;
            do
            {
                X[0] = A[0]+(B[0]-A[0])*t;
                X[1] = A[1]+(B[1]-A[1])*t;
                X[2] = A[2]+(B[2]-A[2])*t;

                phiX = phi.eval(X);
                phi.eval_grad(X, grad_phi);

                dphidt = grad_phi[0]*(B[0]-A[0])+
                         grad_phi[1]*(B[1]-A[1])+
                         grad_phi[2]*(B[2]-A[2]);
                
                t -= phiX/dphidt;

                it += 1;
            }
            while ((it < n_it) && (std::abs(phiX) > this->tol));

            if (std::abs(phiX) > this->tol)
            {
                io::warning("mesh.hpp - split_edges",
                            "Could not find the intersection between the distance function and the edge.");
            }

            // ADD A NEW VERTEX
            VoronoiVertex new_vertex;
        }
        // ------------------------------------------------------------
    }
    
    void split_faces();
    void split_cells();
    // ================================================================

    // ADD WALL =======================================================
    template <class WALL>
    void add_wall(const WALL & wall)
    {
        // PARAMETERS -------------------------
        const int n_cells = this->cells.size();
        // ------------------------------------

        // VARIABLES --------------------------------------------------
        std::vector<int> cells_to_be_removed, cells_to_be_cut;
        // ------------------------------------------------------------

        // LOOP OVER THE CELLS AND DETERMINE WHICH CELLS WILL BE: -----
        // 1) REMOVED
        // 2) CUT
        // 3) LEFT UNCHANGED
        // ------------------------------------------------------------
        for (int c = 0; c < n_cells; ++c)
        {
            if (this->vc_ptrs[c] != nullptr)
            {
                const VoronoiCell & cell = this->cells[c];
                const int n_vertices = cell.v_conn.size();

                // LOCATION OF THE FIRST VERTEX WITH RESPECT TO THE
                // WALL
                const double init_distance = wall.eval(this->vertices[cell.v_conn[0]].x);

                // DETERMINE WHETHER THERE IS A CHANGE OF SIGN WHILE
                // LOOPING OVER THE VERTICES
                bool same_sign = true;
                int v = 1;
                while (same_sign && (v < n_vertices))
                {
                    const double vertex_distance = wall.eval(this->vertices[cell.v_conn[v]].x);

                    if (init_distance*vertex_distance < 0.0)
                    {
                        same_sign = false;
                    }
                    else
                    {
                        v += 1;
                    }
                }

                // CELL TO BE REMOVED
                if (same_sign && (init_distance > 0.0))
                {
                    cells_to_be_removed.push_back(c);
                }
                // CELL TO BE LEFT UNCHANGED
                else if (same_sign && (init_distance < 0.0))
                {
                    // Do nothing
                }
                // CELL TO BE CUT
                else
                {
                    cells_to_be_cut.push_back(c);
                }
            }
        }
        // ------------------------------------------------------------

        // HANDLE THE DISAPPEARING OF THE CELLS TO BE REMOVED ---------
        for (size_t ic = 0; ic < cells_to_be_removed.size(); ++ic)
        {
            const int c = cells_to_be_removed[ic];

            // TEMPORARILY UPDATE THE NEIGHBOR INFORMATION
            std::vector<int> nbr;
            this->vc_ptrs[c]->neighbors(nbr);
            const int n_nbr = nbr.size();
            for (int n = 0; n < n_nbr; ++n)
            {
                const int nbr_id = nbr[n];

                // Skip the walls
                if (nbr_id < 0)
                {
                    continue;
                }

                // Update remaining
                const int nbr_c = this->get_cell_index(nbr_id);
                voro::voronoicell_neighbor * nbr_vc = this->vc_ptrs[nbr_c];

                for (int v = 0; v < nbr_vc->p; ++v)
                {
                    for (int f = 0; f < nbr_vc->nu[v]; ++f)
                    {
                        if (nbr_vc->ne[v][f] == this->cells[c].id)
                        {
                            nbr_vc->ne[v][f] = -10001;
                        }
                    }  
                }
                nbr_vc->check_relations();
            }

            // REMOVE THE CELLS
            delete this->vc_ptrs[c];
            this->vc_ptrs[c] = nullptr;
        }
        // ------------------------------------------------------------

        // CUT THE CELLS TO BE CUT ------------------------------------
        for (size_t ic = 0; ic < cells_to_be_cut.size(); ++ic)
        {
            const int c = cells_to_be_cut[ic];
            const double * S = &this->seeds[3*c];
            double dist, un[3];
            dist = wall.eval(S);
            wall.eval_grad(S, un);
            const double tmp = 1.0/std::sqrt(un[0]*un[0]+un[1]*un[1]+un[2]*un[2]);
            un[0] *= tmp;
            un[1] *= tmp;
            un[2] *= tmp;

            this->vc_ptrs[c]->nplane(un[0], un[1], un[2], -2.0*dist, wall.id);
            this->vc_ptrs[c]->check_relations();
        }
        // ------------------------------------------------------------

        // UPDATE THE TESSELLATION ------------------------------------
        std::vector<int> ids(n_cells);
        for (int i = 0; i < n_cells; ++i)
        {
            ids[i] = this->cells[i].id;
        }
        this->init_voronoi_cells(ids);
        this->init_voronoi_faces();
        this->init_voronoi_vertices_and_edges();
        this->init_voronoi_walls();
        // ------------------------------------------------------------
    }
    // ================================================================

    // ADD CRACK ======================================================
    template <class CRACK>
    void add_crack(const CRACK & crack)
    {
        // PARAMETERS -------------------------
        const int n_cells = this->cells.size();
        // ------------------------------------

        // VARIABLES -----------------------------------------
        std::vector<int> cells_to_be_split;
        std::vector<voro::voronoicell_neighbor *> new_vc_ptrs;
        // ---------------------------------------------------

        // LOOP OVER THE CELLS AND DETERMINE WHICH CELLS WILL BE SPLIT
        for (int c = 0; c < n_cells; ++c)
        {
            if (this->vc_ptrs[c] != nullptr)
            {
                const VoronoiCell & cell = this->cells[c];
                const int n_vertices = cell.v_conn.size();

                // LOCATION OF THE FIRST VERTEX WITH RESPECT TO THE
                // SURFACE DEFINING THE CRACK
                const double init_distance = crack.eval(this->vertices[cell.v_conn[0]].x);

                // DETERMINE WHETHER THERE IS A CHANGE OF SIGN WHILE
                // LOOPING OVER THE VERTICES
                bool same_sign = true;
                int v = 1;
                while (same_sign && (v < n_vertices))
                {
                    const double vertex_distance = crack.eval(this->vertices[cell.v_conn[v]].x);

                    if (init_distance*vertex_distance < 0.0)
                    {
                        same_sign = false;
                    }
                    else
                    {
                        v += 1;
                    }
                }

                // LOCATION OF THE FIRST VERTEX WITH RESPECT TO THE
                // CRACK TIP
                const double init_distance_from_crack_tip = crack.eval_crack_tip(this->vertices[cell.v_conn[0]].x);

                // DETERMINE WHETHER THE SIGN REMAINS NEGATIVE WHILE
                // LOOPING OVER THE VERTICES
                bool negative_sign = (init_distance_from_crack_tip < 0.0);
                v = 1;
                while (negative_sign && (v < n_vertices))
                {
                    const double vertex_distance = crack.eval_crack_tip(this->vertices[cell.v_conn[v]].x);

                    if (vertex_distance > 0.0)
                    {
                        negative_sign = false;
                    }
                    else
                    {
                        v += 1;
                    }
                }

                if (!same_sign && negative_sign)
                {
                    cells_to_be_split.push_back(c);
                }
            }
        }

        const int n_new_cells = cells_to_be_split.size();
        // ------------------------------------------------------------

        // ADD A COPY OF THE CELLS TO BE SPLIT AND STORE THEIR IDS ----
        std::vector<int> old_ids(n_new_cells, 0);
        for (int ic = 0; ic < n_new_cells; ++ic)
        {
            const int c = cells_to_be_split[ic];
            const VoronoiCell & split_cell = this->cells[c];

            old_ids[ic] = split_cell.id;

            new_vc_ptrs.push_back(new voro::voronoicell_neighbor());
            *new_vc_ptrs[ic] = *this->vc_ptrs[c];
        }
        // ------------------------------------------------------------

        // EVAL A SET OF NEW IDS FOR THE NEW CELLS --------------------
        std::vector<int> new_ids(n_new_cells, 0);
        if (n_new_cells > 0)
        {
            new_ids[0] = this->get_new_id();
            for (int ic = 1; ic < n_new_cells; ++ic)
            {
                new_ids[ic] = new_ids[ic-1]+1;
            }
        }
        // ------------------------------------------------------------

        // SPLIT THE CELLS TO BE SPLIT --------------------------------
        for (int ic = 0; ic < n_new_cells; ++ic)
        {
            const int c = cells_to_be_split[ic];

            const double * S = &this->seeds[3*c];
            const double dist = crack.eval(S);
            double un[3];
            crack.eval_grad(S, un);
            const double tmp = 1.0/std::sqrt(un[0]*un[0]+un[1]*un[1]+un[2]*un[2]);
            un[0] *= tmp;
            un[1] *= tmp;
            un[2] *= tmp;

            this->vc_ptrs[c]->nplane(un[0], un[1], un[2], -2.0*dist, new_ids[ic]);
            this->vc_ptrs[c]->check_relations();

            new_vc_ptrs[ic]->nplane(-un[0], -un[1], -un[2], +2.0*dist, old_ids[ic]);
            new_vc_ptrs[ic]->check_relations();
        }
        // ------------------------------------------------------------

        // UPDATE NEIGHBORS' INFO -------------------------------------
        for (int ic = 0; ic < n_new_cells; ++ic)
        {
            const int c = cells_to_be_split[ic];

            voro::voronoicell_neighbor * cracked_vc = this->vc_ptrs[c];

            // Locate all neighbors touching the cracked cells
            int ii, jj;
            if (voro_utils::find_face(*cracked_vc, ii, jj, new_ids[ic]))
            {
                io::error("mesh.hpp - add_crack",
                          "Error locating cut face");
            }
            std::vector<int> cut_nbr;
            voro_utils::face_neighbors(*cracked_vc, ii, jj, new_ids[ic], cut_nbr);
            const int n_cut_nbr = cut_nbr.size();

            // Loop over the neighbors and split the affected faces
            for (int cn = 0; cn < n_cut_nbr; ++cn)
            {
                // Skip the walls
                if (cut_nbr[cn] < 0)
                {

                }
                else
                {
                    const int cut_nbr_c = this->get_cell_index(cut_nbr[cn]);
                    std::vector<int>::iterator it;

                    it = std::find(cells_to_be_split.begin(), cells_to_be_split.end(), cut_nbr_c);

                    if (it == cells_to_be_split.end())
                    {
                        voro::voronoicell_neighbor * nbr_vc = this->vc_ptrs[cut_nbr_c];

                        const double * nbr_S = &this->seeds[3*cut_nbr_c];
                        const double nbr_dist = crack.eval(nbr_S);
                        double nbr_un[3];
                        crack.eval_grad(nbr_S, nbr_un);
                        const double tmp = 1.0/std::sqrt(nbr_un[0]*nbr_un[0]+nbr_un[1]*nbr_un[1]+nbr_un[2]*nbr_un[2]);
                        nbr_un[0] *= tmp;
                        nbr_un[1] *= tmp;
                        nbr_un[2] *= tmp;

                        // Perform the splitting and check the edge relation info is
                        // correct
                        voro_utils::split_face(*nbr_vc, old_ids[ic], new_ids[ic],
                                               nbr_un[0], nbr_un[1], nbr_un[2], -2.0*nbr_dist);
                        nbr_vc->check_relations();
                    }
                }
            }

            // Update the neighboring information of the faces of the
            // other voro++ cell
            {
                std::vector<int>::iterator it;

                for (int v = 0; v < new_vc_ptrs[ic]->p; ++v)
                {
                    for (int f = 0; f < new_vc_ptrs[ic]->nu[v]; ++f)
                    {
                        if (new_vc_ptrs[ic]->ne[v][f] == old_ids[ic])
                        {
                            continue;
                        }

                        it = std::find(old_ids.begin(), old_ids.end(), new_vc_ptrs[ic]->ne[v][f]);
                        
                        if (it != old_ids.end())
                        {
                            new_vc_ptrs[ic]->ne[v][f] = new_ids[std::distance(old_ids.begin(), it)];
                        }
                    }  
                }
                new_vc_ptrs[ic]->check_relations();
            }

            // Get the neighbors of the other voro++ cell
            std::vector<int> nbr;
            new_vc_ptrs[ic]->neighbors(nbr);

            // Loop over the neighbors of the other voro++ cell and
            // update the neighboring information of the faces
            const int n_nbr = nbr.size();
            for (int n = 0; n < n_nbr; ++n)
            {
                const int nbr_id = nbr[n];

                std::vector<int>::iterator it;

                // Skip the walls or the original cell
                if ((nbr_id < 0) || (nbr_id == old_ids[ic]))
                {
                    continue;
                }

                // Skip those cells whose faces have been cut
                bool found = false;
                for (int cn = 0; cn < n_cut_nbr; ++cn)
                {
                    if (nbr_id == cut_nbr[cn])
                    {
                        found = true;
                    }
                }
                if (found)
                {
                    continue;
                }

                // Skip the newly added cells
                it = std::find(new_ids.begin(), new_ids.end(), nbr_id);
                
                if (it != new_ids.end())
                {
                    continue;
                }

                // Update remaining
                const int nbr_c = this->get_cell_index(nbr_id);
                voro::voronoicell_neighbor * nbr_vc = this->vc_ptrs[nbr_c];

                for (int v = 0; v < nbr_vc->p; ++v)
                {
                    for (int f = 0; f < nbr_vc->nu[v]; ++f)
                    {   
                        if (nbr_vc->ne[v][f] == old_ids[ic])
                        {
                            nbr_vc->ne[v][f] = new_ids[ic];
                        }
                    }  
                }
                nbr_vc->check_relations();
            }
        }
        // ------------------------------------------------------------

        // UPDATE THE vc_ptr VECTOR -----------------------------------
        for (int ic = 0; ic < n_new_cells; ++ic)
        {
            this->vc_ptrs.push_back(new_vc_ptrs[ic]);
        }
        // ------------------------------------------------------------

        // UPDATE THE seeds VECTOR ------------------------------------
        for (int ic = 0; ic < n_new_cells; ++ic)
        {
            const int c = cells_to_be_split[ic];
            const double * S = &this->seeds[3*c];

            this->seeds.push_back(S[0]);
            this->seeds.push_back(S[1]);
            this->seeds.push_back(S[2]);
        }
        // ------------------------------------------------------------

        // UPDATE TOLERANCE -------------------------------------------
        this->init_tolerance();
        // ------------------------------------------------------------

        // RECOMPUTE THE VORONOI TESSELLATION -------------------------
        {
            std::vector<int> ids(n_cells+n_new_cells);
            for (int i = 0; i < n_cells; ++i)
            {
                ids[i] = this->cells[i].id;
            }
            for (int ic = 0; ic < n_new_cells; ++ic)
            {
                ids[n_cells+ic] = new_ids[ic];
            }
            this->init_voronoi_cells(ids);
            this->init_voronoi_faces();
            this->init_voronoi_vertices_and_edges();
            this->init_voronoi_walls();
        }
        // ------------------------------------------------------------

        // UPDATE THE CRACKS INFO -------------------------------------
        for (int ic = 0; ic < n_new_cells; ++ic)
        {
            const VoronoiCell & cell = this->cells[n_cells+ic];
            const int n_cell_faces = cell.f_conn.size();

            for (int cf = 0; cf < n_cell_faces; ++cf)
            {
                const int f = cell.f_conn[cf];
                VoronoiFace & face = this->faces[f];

                if (cell.nbr[cf] == cells_to_be_split[ic])
                {
                    face.bou_type = BOU_TYPE_CRACK;
                }
            }
        }
        // ------------------------------------------------------------
    }
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