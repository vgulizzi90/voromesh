// main.cpp

// INCLUDES ===========================================================
#include <iostream>
#include <voro++.hh>
#include <voromesh++.hpp>
// ====================================================================



// RANDOM NUMBER IN THE INTERVAL [0,1] ================================
double rnd()
{
    return (double(rand())/RAND_MAX);
}
// ====================================================================


// MAIN ===============================================================
int main()
{
    // CONSTANTS ------------------------------------------------------
    // Container geometry parameters
    const double X1L = -1.0;
    const double X1U = +1.0;
    const double X2L = -1.0;
    const double X2U = +1.0;
    const double X3L = -1.0;
    const double X3U = +1.0;

    // Container blocks
    const int nb1 = 6;
    const int nb2 = 6;
    const int nb3 = 6;

    // Periodic flags
    const bool p1 = false;
    const bool p2 = false;
    const bool p3 = false;

    // Number of particles
    const int n_particles = 1000;
    // ----------------------------------------------------------------

    // VARIABLES ------------------------------------------------------
    double X1, X2, X3;

    voro::container con(X1L, X1U, X2L, X2U, X3L, X3U,
                        nb1, nb2, nb3,
                        p1, p2, p3,
                        8);
    
    voromesh::Mesh msh;
    // ----------------------------------------------------------------

    // PLACE THE PARTICLES --------------------------------------------
    for (int id = 0; id < n_particles; ++id)
    {
        X1 = X1L+rnd()*(X1U-X1L);
        X2 = X2L+rnd()*(X2U-X2L);
        X3 = X3L+rnd()*(X3U-X3L);
        con.put(id, X1, X2, X3);
    }
    // ----------------------------------------------------------------

    // INITIALIZE THE MESH WITH THE TESSELLATION INFORMATION
    msh.init(con);

    // EXPORT TO VTK FORMAT
    msh.export_tessellation_vtk("tess.vtu");
    // -----------------------------------------------------

    // REPORT ---------------------------------------------------------
    
    // ----------------------------------------------------------------

    // FILL IN THE OUTPUT DATA STRUCTURES -----------------------------
    std::vector<int> cid(n_particles);
    std::vector<std::array<double, 4>> cc(n_particles);
    std::vector<std::vector<double>> cn(n_particles);
    
    for (int i = 0; i < n_particles; ++i)
    {
        cid[i] = msh.cells[i].id;

        cc[i][0] = msh.cells[i].centroid[0];
        cc[i][1] = msh.cells[i].centroid[1];
        cc[i][2] = msh.cells[i].centroid[2];
        cc[i][3] = msh.cells[i].volume;

        const int n_nbr = msh.cells[i].nbr.size();
        cn[i].push_back(n_nbr);

        for (int n = 0; n < n_nbr; ++n)
        {
            const int nbr_c = msh.cells[i].nbr[n];
        }
    }
    // ----------------------------------------------------------------
}
// ====================================================================