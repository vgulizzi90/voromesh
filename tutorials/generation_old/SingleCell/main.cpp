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



// MAIN ===============================================================
int main()
{
    // CONSTANTS ------------------------------------------------------
    // ID of the particles
    const int id = 0;

    // Number of planes
    const int n_planes = 250;
    // ----------------------------------------------------------------

    // VARIABLES ------------------------------------------------------
    double X1, X2, X3, r, rsq;
    
    voro::voronoicell_neighbor vc;
    
    voromesh::Mesh msh;
    // ----------------------------------------------------------------

    // INITIALIZE THE VORO++ CELL AS A CUBE ---------------------------
    vc.init(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
    // ----------------------------------------------------------------

    // CUT THE CELL BY n_planes RANDOM PLANES -------------------------
    for (int i = 0; i < n_planes; ++i)
    {
        X1 = 2.0*rnd()-1.0;
        X2 = 2.0*rnd()-1.0;
        X3 = 2.0*rnd()-1.0;
        rsq = X1*X1+X2*X2+X3*X3;
        if ((rsq > 0.01) && (rsq < 1.0))
        {
            r = 1.0/std::sqrt(rsq);
            X1 *= r;
            X2 *= r;
            X3 *= r;
            vc.nplane(X1, X2, X3, 1.0, -(i+1));
        }
    }
    // ----------------------------------------------------------------

    // ADD THE CELL ------------
    msh.add_vc(vc, id);

    // BUILD THE MESH
    msh.build(3, 0.05);

    // EXPORT TO VTK FORMAT
    msh.export_vtk("mesh.vtu");
    // ------------------------

    // PERFORM FIRST CUT --------------------------
    {
        const double un[3] = {rnd(), rnd(), rnd()};
        msh.cut_cell_by_vector(id, un);
    }

    // BUILD THE MESH
    msh.build(3, 0.05);

    // EXPORT TO VTK FORMAT
    msh.export_vtk("mesh-2.vtu");
    // --------------------------------------------

    // PERFORM SECOND CUT ON THE NEWLY ADDED CELL -
    {
        const double un[3] = {rnd(), rnd(), rnd()};
        msh.cut_cell_by_vector(id+1, un);
    }

    // BUILD THE MESH
    msh.build(3, 0.05);

    // EXPORT TO VTK FORMAT
    msh.export_vtk("mesh-3.vtu");
    // --------------------------------------------
}
// ====================================================================