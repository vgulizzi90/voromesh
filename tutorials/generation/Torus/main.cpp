// main.cpp

// INCLUDES ===========================================================
#include <iostream>
#include <voro++.hh>
#include <voromesh++.hpp>

#include <cmath>
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

    // CHECK WHETHER A POINT IS INSIDE ================================
    bool point_inside(double X1, double X2, double X3)
    {
        const double RX = std::sqrt(X1*X1+X2*X2);
        const double r2 = this->r*this->r;
        const double tmp = RX-this->R;
        return ((tmp*tmp+X3*X3) < r2);
    }
    // ================================================================

    // WALL CUT =======================================================
    template <class VC>
    inline bool cut_cell_base(VC & vc, double X1, double X2, double X3)
    {
        const double RX = std::sqrt(X1*X1+X2*X2);
        double ds = RX-this->R;
        double dst = ds*ds+X3*X3;
        if (dst > 0.01*this->r)
        {
            dst = 2.0*this->r/std::sqrt(dst)-2.0;
            X3 *= dst;
            ds *= dst/RX;
            X1 *= ds;
            X2 *= ds;
            return vc.nplane(X1, X2, X3, this->id);
        }
        return true;
    }

    bool cut_cell(voro::voronoicell & vc, double X1, double X2, double X3)
    {
        return cut_cell_base(vc, X1, X2, X3);
    }
    bool cut_cell(voro::voronoicell_neighbor & vc, double X1, double X2, double X3)
    {
        return cut_cell_base(vc, X1, X2, X3);
    }
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
    const double L = oR+0.5;
    const double H = r+0.5;

    const double X1L = -L;
    const double X1U = +L;
    const double X2L = -L;
    const double X2U = +L;
    const double X3L = -H;
    const double X3U = +H;

    // Container blocks
    const int nb1 = 10;
    const int nb2 = 10;
    const int nb3 = 3;

    // Periodic flags
    const bool p1 = false;
    const bool p2 = false;
    const bool p3 = false;
    // ----------------------------------------------------------------

    // VARIABLES ------------------------------------------------------
    voro::container con(X1L, X1U, X2L, X2U, X3L, X3U,
                        nb1, nb2, nb3,
                        p1, p2, p3,
                        8);
    
    voromesh::Mesh msh;
    // ----------------------------------------------------------------

    // ADD THE WALL TO THE CONTAINER
    wall_torus tor(R, r);
    con.add_wall(tor);
    // -----------------------------

    // IMPORT PARTICLES FROM FILE
	con.import("pack_torus.txt");
    // --------------------------

    // INITIALIZE THE MESH WITH THE TESSELLATION INFORMATION
    msh.init(con);

    // BUILD THE MESH
    msh.build(3, 0.4);

    // EXPORT TO VTK FORMAT
    msh.export_vtk("mesh.vtu");
    // -----------------------------------------------------
}
// ====================================================================