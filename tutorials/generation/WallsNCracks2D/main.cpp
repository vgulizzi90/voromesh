// main.cpp

// INCLUDES ===========================================================
#include <algorithm>

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



// WALL CLASS #########################################################
struct plane_wall
{
    // DATA MEMBERS ===================================================
    const double un[3];
    const double d;
    const int id;
    // ================================================================

    // CONSTRUCTOR ====================================================
    plane_wall(const double * un_, const double d_, const int id_ = -99)
    :
    un{un_[0], un_[1], un_[2]},
    d(d_),
    id(id_)
    {}
    // ================================================================

    // EVAL THE DISTANCE FUNCTION DEFINING THE WALL ===================
    double eval(const double * X) const
    {
        const double res = un[0]*X[0]+un[1]*X[1]+un[2]*X[2]+d;
        return res;
    }

    void eval_grad(const double * X, double * grad) const
    {
        grad[0] = un[0];
        grad[1] = un[1];
        grad[2] = un[2];
    }
    // ================================================================
};
// ####################################################################



// CRACK CLASS ########################################################
struct plane_elliptic_crack
{
    // DATA MEMBERS ===================================================
    const double un[3];
    const double d;
    const double C[3], R[3];
    const int id;
    // ================================================================

    // CONSTRUCTOR ====================================================+
    plane_elliptic_crack(const double * un_, const double d_,
                         const double * C_, const double * R_,
                         const int id_ = -1099)
    :
    un{un_[0], un_[1], un_[2]},
    d(d_),
    C{C_[0], C_[1], C_[2]},
    R{R_[0], R_[1], R_[2]},
    id(id_)
    {}
    // ================================================================

    // EVAL THE DISTANCE FUNCTIONS DEFINING THE CRACK =================
    double eval(const double * X) const
    {
        const double res = un[0]*X[0]+un[1]*X[1]+un[2]*X[2]+d;
        return res;
    }

    void eval_grad(const double * X, double * grad) const
    {
        grad[0] = un[0];
        grad[1] = un[1];
        grad[2] = un[2];
    }

    double eval_crack_tip(const double * X) const
    {
        const double x[3] = {(X[0]-C[0])/R[0], (X[1]-C[1])/R[1], (X[2]-C[2])/R[2]};
        const double res = x[0]*x[0]+x[1]*x[1]+x[2]*x[2]-1.0;
        return res;
    }
    
    void eval_grad_crack_tip(const double * X, double * grad) const
    {
        grad[0] = 2.0*(X[0]-C[0])/(R[0]*R[0]);
        grad[1] = 2.0*(X[1]-C[1])/(R[1]*R[1]);
        grad[2] = 2.0*(X[2]-C[2])/(R[2]*R[2]);
    }
    // ================================================================
};
// ####################################################################



// MAIN ===============================================================
int main()
{
    // CONSTANTS ------------------------------------------------------
    // Container geometry parameters
    const double X1L = -1.0;
    const double X1U = +1.0;
    const double X2L = -1.0;
    const double X2U = +1.0;
    const double X3L = -0.025;
    const double X3U = +0.025;

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

    // PLACE THE PARTICLES (TO GENERATE A COLUMNAR MORPHOLOGY) --------
    for (int i = 0; i < n_particles; ++i)
    {
        X1 = X1L+rnd()*(X1U-X1L);
        X2 = X2L+rnd()*(X2U-X2L);
        X3 = 0.5*(X3U+X3L);
        con.put(i, X1, X2, X3);
    }
    // ----------------------------------------------------------------

    // INITIALIZE THE MESH WITH THE TESSELLATION INFORMATION
    msh.init(con);

    // BUILD THE MESH
    msh.build(3, 0.1);

    // EXPORT TO VTK FORMAT
    msh.export_vtk("mesh.vtu");
    // -----------------------------------------------------

    // ADD WALLS ------------------------------------------------------
    {
        const double un[3] = {0.0, 1.0, 0.0};
        const double d = -0.6;
        const plane_wall pln(un, d, -101);
        msh.add_wall(pln);

        // BUILD THE MESH
        msh.build(3, 0.1);

        // EXPORT TO VTK FORMAT
        msh.export_vtk("mesh-2.vtu");
    }

    {
        const double alpha = M_PI/3.0;
        const double un[3] = {std::cos(alpha), std::sin(alpha), 0.0};
        const double d = -std::sin(alpha)*0.6;
        const plane_wall pln(un, d, -102);
        msh.add_wall(pln);

        // BUILD THE MESH
        msh.build(3, 0.1);

        // EXPORT TO VTK FORMAT
        msh.export_vtk("mesh-3.vtu");
    }
    // ----------------------------------------------------------------

    // ADD CRACKS -----------------------------------------------------
    {
        const double alpha = M_PI/2.0;
        const double un[3] = {std::cos(alpha), std::sin(alpha), 0.0};
        const double C[3] = {-0.5, -0.2, 0.0};
        const double R[3] = {0.4, 0.4, 0.4};
        const double d = -(C[0]*un[0]+C[1]*un[1]+C[2]*un[2]);
        const plane_elliptic_crack crk(un, d, C, R, -1002);
        msh.add_crack(crk);

        // BUILD THE MESH
        msh.build(3, 0.1);

        // EXPORT TO VTK FORMAT
        msh.export_vtk("mesh-4.vtu");
    }

    {
        const double alpha = M_PI/3.0;
        const double un[3] = {std::cos(alpha), std::sin(alpha), 0.0};
        const double C[3] = {-0.5, -0.2, 0.0};
        const double R[3] = {0.4, 0.4, 0.4};
        const double d = -(C[0]*un[0]+C[1]*un[1]+C[2]*un[2]);
        const plane_elliptic_crack crk(un, d, C, R, -1003);
        msh.add_crack(crk);

        // BUILD THE MESH
        msh.build(3, 0.1);

        // EXPORT TO VTK FORMAT
        msh.export_vtk("mesh-5.vtu");
    }
    // ----------------------------------------------------------------

    // REPORT WALLS ---------------------------------------------------
    // ----------------------------------------------------------------

    // REPORT CRACKS --------------------------------------------------
    // ----------------------------------------------------------------
}
// ====================================================================