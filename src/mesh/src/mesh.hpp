// mesh.hpp

#ifndef VOROMESH_MESH_HPP_
#define VOROMESH_MESH_HPP_

#include <vector>
#include <algorithm>

#include <voro++.hh>

#include "utils_io.hpp"
#include "utils_math.hpp"
#include "utils_voropp.hpp"

#define ID_UNDEFINED -1

#define BOU_TYPE_UNDEFINED -1
#define BOU_TYPE_WALL 0
#define BOU_TYPE_INTERFACE 1
#define BOU_TYPE_PERIODIC 2
#define BOU_TYPE_CRACK 3

#define BOU_TYPE_BAD_NEIGHBOR -1010
#define WALL_ID_BAD_NEIGHBOR -1010

namespace voromesh
{
// VORONOI FACE CLASS #################################################
struct VoronoiFace
{
    // DATA MEMBERS ===================================================
    int bou_type;
    int parent_cells[2];
    double centroid[3], area;
    // ================================================================

    // CONSTRUCTOR ====================================================
    VoronoiFace();
    // ================================================================

    // EVAL LOCAL REFERENCE SYSTEM ====================================
    // ================================================================

    // BUILD MESH =====================================================
    // ================================================================
};
// ####################################################################



// VORONOI CELL CLASS #################################################
struct VoronoiCell
{
    // DATA MEMBERS ===================================================
    int id;
    std::vector<int> nbr;
    double centroid[3], volume;

    std::vector<int> i_conn;

    std::vector<int> f_conn;
    // ================================================================

    // CONSTRUCTOR ====================================================
    VoronoiCell();
    // ================================================================

    // BUILD MESH =====================================================
    // ================================================================
};
// ####################################################################



// VORONOI WALL CLASS #################################################
struct VoronoiWall
{
    // DATA MEMBERS ===================================================
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
    int id;

    std::vector<int> f_conn;
    std::vector<std::array<int, 2>> c_conn;
    // ================================================================

    // CONSTRUCTOR ====================================================
    VoronoiCrack();
    // ================================================================

    // AUXILIARY METHODS ==============================================
    bool contains_cell(const int c) const;
    // ================================================================
};
// ####################################################################


    
// MESH CLASS #########################################################
struct Mesh
{
    // DATA MEMBERS ===================================================
    std::vector<double> seeds;
    std::vector<int> ids;

    std::vector<VoronoiCell> cells;
    std::vector<VoronoiFace> faces;
    std::vector<VoronoiWall> walls;
    std::vector<VoronoiCrack> cracks;

    // VORO++
    std::vector<voro::voronoicell_neighbor *> vc_ptrs;

    // ADDITIONAL FLAGS
    bool columnar[3];
    // ================================================================

    // CONSTRUCTOR ====================================================
    Mesh();
    // ================================================================

    // DESTRUCTOR =====================================================
    ~Mesh();
    // ================================================================

    // FREE MEMORY ====================================================
    void clear_tessellation();
    void free();
    // ================================================================

    // INITIALIZATION FROM A VORO++ CONTAINER =========================
    template <typename CONTAINER>
    void init(CONTAINER & con)
    {
        // PARAMETERS -----------------------
        const int np = con.total_particles();
        // ----------------------------------

        // VARIABLES --------------------------------------------------
        double S1_min, S1_max, S2_min, S2_max, S3_min, S3_max;
        double L_max;

        S1_min = std::numeric_limits<double>::max();
        S1_max = std::numeric_limits<double>::min();
        S2_min = std::numeric_limits<double>::max();
        S2_max = std::numeric_limits<double>::min();
        S3_min = std::numeric_limits<double>::max();
        S3_max = std::numeric_limits<double>::min();
        // ------------------------------------------------------------

        // FREE MEMORY
        this->free();
        // -----------
        
        // INIT MEMORY ----------
        this->seeds.resize(3*np);
        this->ids.resize(np);
        this->vc_ptrs.resize(np);
        // ----------------------

        // COMPUTE THE VORONOI CELLS ----------------------------------
        int k = 0;
        double * pp;
        voro::voronoicell_neighbor * vc;
        voro::c_loop_all cl(con);

        if (cl.start()) do
        {
            vc = new voro::voronoicell_neighbor();
            
            if (con.compute_cell(*vc, cl))
            {
                pp = con.p[cl.ijk]+3*cl.q;
                this->seeds[3*k+0] = pp[0];
                this->seeds[3*k+1] = pp[1];
                this->seeds[3*k+2] = pp[2];

                S1_min = std::min(S1_min, this->seeds[3*k+0]);
                S1_max = std::max(S1_max, this->seeds[3*k+0]);
                S2_min = std::min(S2_min, this->seeds[3*k+1]);
                S2_max = std::max(S2_max, this->seeds[3*k+1]);
                S3_min = std::min(S3_min, this->seeds[3*k+2]);
                S3_max = std::max(S3_max, this->seeds[3*k+2]);

                // Check the id of the voronoi cell is not negative
                if (cl.pid() < 0)
                {
                    io::error("mesh.hpp - Mesh::init",
                              "Voronoi cells cannot have negative ids.");
                }
                this->ids[k] = cl.pid();

                this->vc_ptrs[k] = vc;
            }
            else
            {
                delete vc;
                this->ids[k] = ID_UNDEFINED;
                this->vc_ptrs[k] = nullptr;
            }

            k += 1;
        }
        while (cl.inc());
        // ------------------------------------------------------------

        // CHECK WHETHER THE TESSELLATION IS COLUMNAR -----------------
        L_max = std::max(S1_max-S1_min, std::max(S2_max-S2_min, S3_max-S3_min));

        if ((S1_max-S1_min) < L_max*1.0e-12) this->columnar[0] = true;
        if ((S2_max-S2_min) < L_max*1.0e-12) this->columnar[1] = true;
        if ((S3_max-S3_min) < L_max*1.0e-12) this->columnar[2] = true;
        // ------------------------------------------------------------

        // INIT THE VORONOI TESSELLATION DATA STRUCTURES
        this->init_tessellation();
        // ---------------------------------------------
    }
    // ================================================================

    // AUXILIARY METHODS ==============================================
    int get_new_id() const;
    int find_cell_index_by_id(const int id) const;
    int find_wall_index_by_id(const int id) const;
    // ================================================================

    // INITIALIZATION: TESSELLATION ===================================
    void init_tessellation();
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
            const double * S = &this->seeds[3*c];
            voro::voronoicell_neighbor * vc = this->vc_ptrs[c];

            if (vc != nullptr)
            {
                std::vector<double> vertices;
                vc->vertices(S[0], S[1], S[2], vertices);

                const int n_vertices = vertices.size()/3;

                // LOCATION OF THE FIRST VERTEX WITH RESPECT TO THE
                // WALL
                const double V0_distance = wall.eval(vertices.data());

                // DETERMINE WHETHER THERE IS A CHANGE OF SIGN WHILE
                // LOOPING OVER THE VERTICES
                bool same_sign = true;
                int v = 1;
                while (same_sign && (v < n_vertices))
                {
                    const double V_distance = wall.eval(&vertices[3*v]);

                    if (V0_distance*V_distance < 0.0)
                    {
                        same_sign = false;
                    }
                    else
                    {
                        v += 1;
                    }
                }

                // CELL TO BE REMOVED
                if (same_sign && (V0_distance > 0.0))
                {
                    cells_to_be_removed.push_back(c);
                }
                // CELL TO BE LEFT UNCHANGED
                else if (same_sign && (V0_distance < 0.0))
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
                const int nbr_c = this->find_cell_index_by_id(nbr_id);
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

        // RECOMPUTE THE TESSELLATION
        this->init_tessellation();
        // --------------------------
    }
    // ================================================================

    // ADD CRACK ======================================================
    template <class CRACK>
    void add_crack(const CRACK & crack)
    {
        // PARAMETERS -------------------------
        const int n_cells = this->cells.size();
        // ------------------------------------

        // VARIABLES --------------------------------------------------
        std::vector<int> cells_to_be_split;
        std::vector<voro::voronoicell_neighbor *> new_vc_ptrs;

        std::vector<std::array<int, 2>> new_cell_pairs_from_old_cracks;
        // ------------------------------------------------------------

        // LOOP OVER THE CELLS AND DETERMINE WHICH CELLS WILL BE SPLIT
        for (int c = 0; c < n_cells; ++c)
        {
            const double * S = &this->seeds[3*c];
            voro::voronoicell_neighbor * vc = this->vc_ptrs[c];

            if (vc != nullptr)
            {
                std::vector<double> vertices;
                vc->vertices(S[0], S[1], S[2], vertices);

                const int n_vertices = vertices.size()/3;

                // LOCATION OF THE FIRST VERTEX WITH RESPECT TO THE
                // SURFACE DEFINING THE CRACK
                const double V0_distance = crack.eval(vertices.data());

                // DETERMINE WHETHER THERE IS A CHANGE OF SIGN WHILE
                // LOOPING OVER THE VERTICES
                bool same_sign = true;
                int v = 1;
                while (same_sign && (v < n_vertices))
                {
                    const double V_distance = crack.eval(&vertices[3*v]);

                    if (V0_distance*V_distance < 0.0)
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
                const double V0_distance_from_crack_tip = crack.eval_crack_tip(vertices.data());

                // DETERMINE WHETHER THE SIGN REMAINS NEGATIVE WHILE
                // LOOPING OVER THE VERTICES
                bool negative_sign = (V0_distance_from_crack_tip < 0.0);
                v = 1;
                while (negative_sign && (v < n_vertices))
                {
                    const double V_distance = crack.eval_crack_tip(&vertices[3*v]);

                    if (V_distance > 0.0)
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
            const VoronoiCell & cell = this->cells[c];

            voro::voronoicell_neighbor * cracked_vc = this->vc_ptrs[c];

            // Locate all neighbors touching the cracked cell
            int ii, jj;
            if (voro_utils::find_face(*cracked_vc, ii, jj, new_ids[ic]))
            {
                io::error("mesh.hpp - add_crack",
                          "Error locating cut face");
            }
            std::vector<int> cut_nbr;
            voro_utils::face_neighbors(*cracked_vc, ii, jj, new_ids[ic], cut_nbr);
            const int cut_n_nbr = cut_nbr.size();

            // Loop over the neighbors and split the affected faces
            for (int cn = 0; cn < cut_n_nbr; ++cn)
            {
                // Skip the walls
                if (cut_nbr[cn] < 0)
                {

                }
                else
                {
                    const int cut_nbr_c = this->find_cell_index_by_id(cut_nbr[cn]);
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

                        /*
                        // Old and new cell pair
                        const std::array<int, 2> cell_pair = {c, cut_nbr_c};
                        const std::array<int, 2> new_cell_pair = {new_ids[ic], cut_nbr_c};

                        // Check whether we are splitting a crack face
                        const int n_cell_faces = cell.f_conn.size();    
                        int cf = 0;
                        bool is_crack = false;
                        while ((!is_crack) && (cf < n_cell_faces))
                        {
                            const VoronoiFace & face = this->faces[cell.f_conn[cf]];
                            if ((face.parent_cells[0] == cell_pair[0]) &&
                                (face.parent_cells[1] == cell_pair[1]) &&
                                (face.bou_type == BOU_TYPE_CRACK))
                            {
                                is_crack = true;
                            }
                            else
                            {
                                cf += 1;
                            }
                        }

                        // If the face that we are splitting is a crack
                        // face, we must find the parent VoronoiCrack
                        // and introduce a new crack face
                        if (is_crack)
                        {
                            const int n_old_cracks = this->cracks.size();
                            int cr = 0;
                            bool found = false;
                            while ((!found) && (cr < n_old_cracks))
                            {
                                const VoronoiCrack & voro_crack = this->cracks[cr];

                                std::vector<int>::const_iterator f_it = std::find(voro_crack.f_conn.begin(),
                                                                                  voro_crack.f_conn.end(),
                                                                                  cell.f_conn[cf]);
                                if (f_it != voro_crack.f_conn.end())
                                {
                                    found = true;
                                }
                                else
                                {
                                    cr += 1;
                                }
                            }

                            if (!found)
                            {
                                std::string msg = "The information of the crack's faces is not consistent:\n";
                                msg += "| current face: "+std::to_string(cell.f_conn[cf])+";\n";
                                msg += "| old cell pair: "+std::to_string(cell_pair[0])+","+std::to_string(cell_pair[1])+"\n";
                                msg += "| new cell pair: "+std::to_string(new_cell_pair[0])+","+std::to_string(new_cell_pair[1])+"\n";
                                io::error("mesh.hpp - Mesh::add_crack", msg);
                            }
                            
                            this->cracks[cr].f_conn.push_back(cell.f_conn[cf]);
                            this->cracks[cr].c_conn.push_back(new_cell_pair);
                        }
                        */
                    }
                }
            }

            // Add crack faces that might have originated from the
            // cutting of old crack faces
            {
                std::vector<int>::iterator it;
                std::vector<int> nbr;
                new_vc_ptrs[ic]->neighbors(nbr);
                const int n_nbr = nbr.size();

                for (int n = 0; n < n_nbr; ++n)
                {
                    if ((nbr[n] < 0) || (nbr[n] == old_ids[ic]))
                    {
                        continue;
                    }

                    const int nbr_c = this->find_cell_index_by_id(nbr[n]);
                    it = std::find(old_ids.begin(), old_ids.end(), nbr[n]);

                    bool do_continue = false;
                    std::array<int, 2> cell_pair;
                    std::array<int, 2> new_cell_pair;

                    if (it != old_ids.end())
                    {
                        if (c > nbr_c)
                        {
                            do_continue = true;
                        }
                        else
                        {
                            cell_pair[0] = c;
                            cell_pair[1] = nbr_c;

                            new_cell_pair[0] = new_ids[ic];
                            new_cell_pair[1] = new_ids[std::distance(old_ids.begin(), it)];
                        }
                    }
                    else
                    {
                        cell_pair[0] = (c < nbr_c) ? c : nbr_c;
                        cell_pair[1] = (c < nbr_c) ? nbr_c : c;

                        new_cell_pair[0] = (c < nbr_c) ? new_ids[ic] : nbr_c;
                        new_cell_pair[1] = (c < nbr_c) ? nbr_c : new_ids[ic];
                    }

                    if (do_continue) continue;

                    {
                        // Check whether we are splitting a crack face
                        const int n_cell_faces = cell.f_conn.size();    
                        int cf = 0;
                        bool is_crack = false;
                        while ((!is_crack) && (cf < n_cell_faces))
                        {
                            const VoronoiFace & face = this->faces[cell.f_conn[cf]];
                            if ((face.parent_cells[0] == cell_pair[0]) &&
                                (face.parent_cells[1] == cell_pair[1]) &&
                                (face.bou_type == BOU_TYPE_CRACK))
                            {
                                is_crack = true;
                            }
                            else
                            {
                                cf += 1;
                            }
                        }

                        // If the face that we are splitting is a crack
                        // face, we must find the parent VoronoiCrack
                        // and introduce a new crack face
                        if (is_crack)
                        {
                            const int n_old_cracks = this->cracks.size();
                            int cr = 0;
                            bool found = false;
                            while ((!found) && (cr < n_old_cracks))
                            {
                                const VoronoiCrack & voro_crack = this->cracks[cr];

                                std::vector<int>::const_iterator f_it = std::find(voro_crack.f_conn.begin(),
                                                                                  voro_crack.f_conn.end(),
                                                                                  cell.f_conn[cf]);
                                if (f_it != voro_crack.f_conn.end())
                                {
                                    found = true;
                                }
                                else
                                {
                                    cr += 1;
                                }
                            }

                            if (!found)
                            {
                                std::string msg = "The information of the crack's faces is not consistent:\n";
                                msg += "| current face: "+std::to_string(cell.f_conn[cf])+";\n";
                                msg += "| old cell pair: "+std::to_string(cell_pair[0])+","+std::to_string(cell_pair[1])+"\n";
                                msg += "| new cell pair: "+std::to_string(new_cell_pair[0])+","+std::to_string(new_cell_pair[1])+"\n";
                                io::error("mesh.hpp - Mesh::add_crack", msg);
                            }
                            
                            this->cracks[cr].f_conn.push_back(cell.f_conn[cf]);
                            this->cracks[cr].c_conn.push_back(new_cell_pair);
                        }
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
                for (int cn = 0; cn < cut_n_nbr; ++cn)
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
                const int nbr_c = this->find_cell_index_by_id(nbr_id);
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

        // RECOMPUTE THE TESSELLATION ---------------------------------
        {
            this->ids.resize(n_cells+n_new_cells);
            for (int c = 0; c < n_cells; ++c)
            {
                this->ids[c] = this->cells[c].id;
            }
            for (int ic = 0; ic < n_new_cells; ++ic)
            {
                this->ids[n_cells+ic] = new_ids[ic];
            }
        }
        this->init_tessellation();
        // ------------------------------------------------------------

        // UPDATE THE OLD FRACTURE DATA STRUCTURES --------------------
        {
            const int n_old_cracks = this->cracks.size();
            
            for (int cr = 0; cr < n_old_cracks; ++cr)
            {
                VoronoiCrack & voro_crack = this->cracks[cr];
                const int n_cell_pairs = voro_crack.c_conn.size();

                int cc = n_cell_pairs-1;

                while (cc >= 0)
                {
                    const int ci = voro_crack.c_conn[cc][0];
                    const int cj = voro_crack.c_conn[cc][1];

                    const VoronoiCell & cell = this->cells[ci];
                    const int n_cell_nbr = cell.nbr.size();

                    int n = 0;
                    bool found = false;

                    while ((!found) && (n < n_cell_nbr))
                    {
                        const int nbr_c = cell.nbr[n];

                        if ((nbr_c >= 0) && (nbr_c == cj))
                        {
                            found = true;
                            const int f = cell.f_conn[n];
                            this->faces[f].bou_type = BOU_TYPE_CRACK;
                            voro_crack.f_conn[cc] = f;
                        }
                        else
                        {
                            n += 1;
                        }
                    }

                    if (!found)
                    {
                        voro_crack.f_conn.erase(voro_crack.f_conn.begin()+cc);
                        voro_crack.c_conn.erase(voro_crack.c_conn.begin()+cc);
                    }

                    cc -= 1;
                }
            }
        }
        // ------------------------------------------------------------

        // ADD A NEW FRACTURE OBJECT ----------------------------------
        {
            VoronoiCrack voro_crack;

            voro_crack.id = crack.id;

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

                        voro_crack.f_conn.push_back(f);
                        voro_crack.c_conn.push_back({face.parent_cells[0], face.parent_cells[1]});
                    }
                }
            }

            this->cracks.push_back(voro_crack);
        }
        // ------------------------------------------------------------
    }
    // ================================================================

    // EXPORT TO VTK FORMAT: TESSELLATION =============================
    void export_tessellation_vtk(const std::string & filepath) const;
    // ================================================================

    // EXPORT TO VTK FORMAT: CENTROIDAL DELAUNAY TRIANGULATION ========
    void export_centroidal_Delaunay(const std::string & filepath) const;
    // ================================================================
};
// ####################################################################

}

#endif