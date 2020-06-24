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

    // Number of fractures
    const int n_fractures = 50;
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

    // EXPORT TO VTK FORMAT
    msh.export_tessellation_vtk("tess.vtu");
    // -----------------------------------------------------

    // ADD FRACTURES --------------------------------------------------
    for (int f = 0; f < n_fractures; ++f)
    {
        const double alpha = 2.0*M_PI*rnd();
        const double un[3] = {std::cos(alpha), std::sin(alpha), 0.0};
        const double C[3] = {X1L+rnd()*(X1U-X1L), X2L+rnd()*(X2U-X2L), 0.5*(X3U+X3L)};
        const double R[3] = {0.4, 0.4, 0.4};
        const double d = -(C[0]*un[0]+C[1]*un[1]+C[2]*un[2]);
        const plane_elliptic_crack crk(un, d, C, R, -(1000+f));

        {
            msh.add_crack(crk);
        }
    }

    msh.export_tessellation_vtk("tess-2.vtu");
    // ----------------------------------------------------------------

    // REPORT WALLS ---------------------------------------------------
    // ----------------------------------------------------------------

    // REPORT CRACKS --------------------------------------------------
    for (int cr = 0; cr < n_fractures; ++cr)
    {
        const voromesh::VoronoiCrack & crack = msh.cracks[cr];
        
        const int n_cell_pairs = crack.c_conn.size();
        
        std::cout << "Crack: " << crack.id << std::endl;
        std::cout << " - cell pairs(" << n_cell_pairs << "):";
        for (int cp = 0; cp < n_cell_pairs; ++cp)
        {
            const int ci = crack.c_conn[cp][0];
            const int cj = crack.c_conn[cp][1];
            std::cout << "(" << msh.cells[ci].id << ", " << msh.cells[cj].id << "), ";
        }
        std::cout << std::endl;
        
    }
    // ----------------------------------------------------------------

    // FILL IN THE OUTPUT DATA STRUCTURES -----------------------------
    // ----------------------------------------------------------------
}
// ====================================================================