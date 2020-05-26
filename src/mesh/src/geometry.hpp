// geometry.hpp

#ifndef VOROMESH_GEOMETRY_HPP_
#define VOROMESH_GEOMETRY_HPP_

#include <vector>
#include <array>
#include <string>
#include <iostream>

#include <voro++.hh>

#include "utils_io.hpp"
#include "utils_math.hpp"

#define BOU_TYPE_WALL 0
#define BOU_TYPE_INTERFACE 1

namespace voromesh
{
    
// VORONOI CELL FACE CLASS ############################################
struct VorocellFace
{
    // DATA MEMBERS ===================================================
    std::vector<double> vertices;
    double centroid[3];
    int parent_cells[2];
    int bou_type;
    // ================================================================

    // DESTRUCTOR =====================================================
    ~VorocellFace();
    // ================================================================
};
// ####################################################################



// VORONOI CELL CLASS #################################################
struct Vorocell
{
    // DATA MEMBERS ===================================================
    voro::voronoicell_neighbor * pp;
    int id;
    double seed[3];
    
    double centroid[3];
    std::vector<int> nbr;
    std::vector<int> f_conn, f_ori;
    // ================================================================

    // DESTRUCTOR =====================================================
    ~Vorocell();
    // ================================================================
};
// ####################################################################



// GEOMETRY CLASS #####################################################
struct Geometry
{
    // DATA MEMBERS ===================================================
    std::vector<Vorocell> cells;
    std::vector<VorocellFace> faces;
    // ================================================================

    // CONSTRUCTOR ====================================================
    Geometry();
    // ================================================================

    // DESTRUCTOR =====================================================
    ~Geometry();
    // ================================================================

    // FREE MEMORY ====================================================
    void free();
    // ================================================================
    
    // INITIALIZATION =================================================
    template <typename CONTAINER>
    void init(CONTAINER & con)
    {
        // PARAMETERS -----------------------
        const int np = con.total_particles();
        // ----------------------------------

        // VARIABLES --------------------
        int k;
        double * pp;
        voro::voronoicell_neighbor * tmp;
        voro::c_loop_all cl(con);
        // ------------------------------

        // FREE MEMORY
        this->free();
        // -----------
        
        // INIT MEMORY
        this->cells.resize(np);

        // COMPUTE THE VORONOI CELLS ----------------------------------
        k = 0;
        if (cl.start()) do
        {
            tmp = new voro::voronoicell_neighbor();
            Vorocell & vc = this->cells[k];
            
            if (con.compute_cell(*tmp, cl))
            {
                pp = con.p[cl.ijk]+3*cl.q;

                vc.pp = tmp;
                vc.id = cl.pid();
                vc.seed[0] = pp[0];
                vc.seed[1] = pp[1];
                vc.seed[2] = pp[2];
            }
            else
            {
                delete tmp;
                vc.pp = nullptr;
                vc.id = -1;
                vc.seed[0] = 0.0;
                vc.seed[1] = 0.0;
                vc.seed[2] = 0.0;
            }

            k += 1;
        }
        while (cl.inc());
        // ------------------------------------------------------------

        // CHECK THAT THE IDS OF THE CELLS ARE NOT NEGATIVE -----------
        for (int k = 0; k < np; ++k)
        {
            Vorocell & vc = this->cells[k];
            if (vc.id < 0)
            {
                io::error("geometry.hpp - Geometry::init",
                          "Voronoi cells cannot have negative ids.");
            }
        }
        // ------------------------------------------------------------

        // BUILD THE GEOMETRICAL ENTITIES
        this->build();
        // ------------------------------
    }
    // ================================================================

    // BUILD THE GEOMETRICAL ENTITIES =================================
    void build();
    // ================================================================

    // EVAL THE GEOMETRIC TOLERANCE ===================================
    double eval_tolerance(const double rel_tol = 1.0e-5);
    // ================================================================

    // EXPORT TO VTK FORMAT ===========================================
    void export_VTK(const std::string & filepath) const;
    // ================================================================
};
// ####################################################################



// GEOMETRY ROUTINES ##################################################

// ####################################################################

}

#endif