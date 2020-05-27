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


// TORUS WALL CLASS ###################################################
struct wall_torus
:
public voro::wall
{
    // DATA MEMBERS ===================================================
    const double R;
    const double r;
    const int id;
    // ================================================================

    // CONSTRUCTOR ====================================================
    wall_torus(const double R_, const double r_, const int id_ = -99)
    :
    R(R_),
    r(r_),
    id(id_)
    {}
    // ================================================================
};
// ####################################################################


// MAIN ===============================================================
int main()
{
    // CONSTANTS ------------------------------------------------------
    // Major and minor torus radii
    const double R = 9.0, r = 3.5;
    
    // Outer radius of the torus
    const double oR = R+r;

    // Container geometry parameters
    const double L = 1.1*oR;
    const double H = 1.1*r;

    const double X1L = -L;
    const double X1U = +L;
    const double X2L = -L;
    const double X2U = +L;
    const double X3L = -H;
    const double X3U = +H;

    // Container blocks
    const int nb1 = 9;
    const int nb2 = 9;
    const int nb3 = 3;

    // Periodic flags
    const bool p1 = false;
    const bool p2 = false;
    const bool p3 = false;

    // Number of particles
    const int n_particles = 734;
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
    msh.build(3, 0.075);

    // EXPORT TO VTK FORMAT
    msh.export_vtk("mesh.vtu");
    // -----------------------------------------------------
}
// ====================================================================