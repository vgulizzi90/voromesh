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



// CRACK CLASS ########################################################
struct elliptic_crack
{
    // DATA MEMBERS ===================================================
    // ================================================================

    // CONSTRUCTOR ====================================================
    // ================================================================

    // EVAL THE DISTANCE FUNCTION =====================================
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
    const int n_particles = 200;
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

    // PERFORM TWO CUTS -----------------------------------------------
    {
        const int id = 128;
        const int c = msh.get_cell_index(id);
        const double S[3] = {msh.seeds[3*c], msh.seeds[3*c+1], msh.seeds[3*c+2]};
        const double un[3] = {rnd(), rnd(), 0.0};
        const double plane[4] = {un[0], un[1], un[2], -(un[0]*S[0]+un[1]*S[1]+un[2]*S[2])};
        msh.cut_cell_by_plane(id, plane);

        const int new_c = n_particles;

        // Find a neighbor that has been affected by the cut
        bool found = false;
        int nbr_c = 0;
        int cnt = 0;
        while ((!found) && (nbr_c < (n_particles+1)))
        {
            const voromesh::VoronoiCell & nbr_cell = msh.cells[nbr_c];
            std::vector<int>::const_iterator it, new_it;

            it = std::find(nbr_cell.nbr.begin(), nbr_cell.nbr.end(), c);
            new_it = std::find(nbr_cell.nbr.begin(), nbr_cell.nbr.end(), new_c);

            if ((it != nbr_cell.nbr.end()) && (new_it != nbr_cell.nbr.end()))
            {
                found = true;
                cnt += 1;
                if (cnt == 1) {found = false; nbr_c += 1;};
            }
            else
            {
                nbr_c += 1;
            }
        }

        // Cut that neighbor as well
        const double new_plane[4] = {un[0], un[1], un[2], 1.0e-12-(un[0]*S[0]+un[1]*S[1]+un[2]*S[2])};
        msh.cut_cell_by_plane(msh.cells[nbr_c].id, new_plane);

        std::cout << "c: " << c << ", id: " << msh.cells[c].id << ", nbr: "; voromesh::math::matrix::print_inline_int(msh.cells[c].nbr.size(), msh.cells[c].nbr.data()); std::cout << std::endl;
        std::cout << "new_c: " << new_c << ", id: " << msh.cells[new_c].id << ", nbr: "; voromesh::math::matrix::print_inline_int(msh.cells[new_c].nbr.size(), msh.cells[new_c].nbr.data()); std::cout << std::endl;
        std::cout << "nbr_c: " << nbr_c << ", id: " << msh.cells[nbr_c].id << ", nbr: "; voromesh::math::matrix::print_inline_int(msh.cells[nbr_c].nbr.size(), msh.cells[nbr_c].nbr.data()); std::cout << std::endl;
    }

    // BUILD THE MESH
    msh.build(3, 0.1);

    // EXPORT TO VTK FORMAT
    msh.export_vtk("mesh-2.vtu");
    // ----------------------------------------------------------------

    // REPORT WALLS ---------------------------------------------------
    {
        const int n_walls = msh.walls.size();
        std::cout << "n_walls: " << n_walls << std::endl;
        for (int w = 0; w < n_walls; ++w)
        {
            const voromesh::VoronoiWall & wall = msh.walls[w];
            const int n_faces = wall.f_conn.size();
            const int n_cells = wall.c_conn.size();
            std::vector<int> cell_ids;

            for (int c = 0; c < n_cells; ++c)
            {
                cell_ids.push_back(msh.cells[wall.c_conn[c]].id);
            }

            std::cout << "wall: " << wall.id << std::endl;
            std::cout << " - faces: "; voromesh::math::matrix::print_inline_int(n_faces, wall.f_conn.data()); std::cout << std::endl;
            std::cout << " - cells: "; voromesh::math::matrix::print_inline_int(n_cells, cell_ids.data()); std::cout << std::endl;
        }
    }
    // ----------------------------------------------------------------

    // REPORT CRACKS --------------------------------------------------
    {
        const int n_cracks = msh.cracks.size();
        std::cout << "n_cracks: " << n_cracks << std::endl;
        for (int cr = 0; cr < n_cracks; ++cr)
        {
            const voromesh::VoronoiCrack & crack = msh.cracks[cr];
            const int n_faces = crack.f_conn.size();
            const int n_cells = crack.c_conn.size();
            std::vector<std::array<int, 2>> cell_pair_ids;

            for (int c = 0; c < n_cells; ++c)
            {
                cell_pair_ids.push_back({msh.cells[crack.c_conn[c][0]].id, msh.cells[crack.c_conn[c][1]].id});
            }

            std::cout << "crack: " << crack.id << std::endl;
            std::cout << " - faces: "; voromesh::math::matrix::print_inline_int(n_faces, crack.f_conn.data()); std::cout << std::endl;
            std::cout << " - cell pairs: "; std::cout << std::endl;
        }
    }
    // ----------------------------------------------------------------
}
// ====================================================================