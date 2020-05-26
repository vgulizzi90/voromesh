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
    const int n_particles = 20;
    // ----------------------------------------------------------------

    // VARIABLES ------------------------------------------------------
    double x1, x2, x3;

    voro::container con(X1L, X1U, X2L, X2U, X3L, X3U,
                        nb1, nb2, nb3,
                        p1, p2, p3,
                        8);
    
    voromesh::Mesh msh;
    // ----------------------------------------------------------------

    // PLACE THE PARTICLES --------------------------------------------
    for (int i = 0; i < n_particles; ++i)
    {
        x1 = X1L+rnd()*(X1U-X1L);
        x2 = X2L+rnd()*(X2U-X2L);
        x3 = X3L+rnd()*(X3U-X3L);
        con.put(i, x1, x2, x3);
    }
    // ----------------------------------------------------------------

    // INITIALIZE THE MESH WITH THE TESSELLATION INFORMATION
    msh.init(con);

    // BUILD THE MESH
    msh.build(3, 1.0);

    // EXPORT TO VTK FORMAT
    msh.export_vtk("mesh.vtu");

    //msh.cut_cell(0, 2.0, 0.0, 0.0, 1.0);
}
// ====================================================================