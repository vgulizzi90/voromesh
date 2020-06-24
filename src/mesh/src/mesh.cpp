// mesh.cpp

#include <limits>
#include "mesh.hpp"
#include "utils_io.hpp"
#include "utils_vtk.hpp"
#include "utils_mesh.hpp"

namespace voromesh
{
// VORONOI FACE CLASS #################################################
// CONSTRUCTOR ========================================================
VoronoiFace::VoronoiFace()
:
bou_type(BOU_TYPE_UNDEFINED),
parent_cells{0, 0},
centroid{0.0, 0.0, 0.0},
area(0.0)
{

}
// ====================================================================
// ####################################################################



// VORONOI CELL CLASS #################################################
// CONSTRUCTOR ========================================================
VoronoiCell::VoronoiCell()
:
id(ID_UNDEFINED),
centroid{0.0, 0.0, 0.0},
volume(0.0)
{

}
// ====================================================================
// ####################################################################



// VORONOI WALL CLASS #################################################
// CONSTRUCTOR ========================================================
VoronoiWall::VoronoiWall()
:
id(ID_UNDEFINED)
{

}
// ====================================================================
// ####################################################################



// VORONOI CRACK CLASS ################################################
// CONSTRUCTOR ========================================================
VoronoiCrack::VoronoiCrack()
:
id(ID_UNDEFINED)
{

}
// ====================================================================

// AUXILIARY METHODS ==================================================
bool VoronoiCrack::contains_cell(const int c) const
{
    const int n_cell_pairs = this->c_conn.size();
    bool found = false;
    int cp = 0;
    while ((!found) && (cp < n_cell_pairs))
    {
        if ((this->c_conn[cp][0] == c) || (this->c_conn[cp][1] == c))
        {
            found = true;
        }
        else
        {
            cp += 1;
        }
    }

    return found;
}
// ====================================================================
// ####################################################################



// MESH CLASS #########################################################
// CONSTRUCTOR ========================================================
Mesh::Mesh()
{

}
// ====================================================================

// DESTRUCTOR =========================================================
Mesh::~Mesh()
{
    this->free();
}
// ====================================================================

// FREE MEMORY ========================================================
void Mesh::clear_tessellation()
{
    this->walls.clear();

    this->faces.clear();
    this->cells.clear();
}

void Mesh::free()
{
    this->clear_tessellation();

    const int n_vcs = this->vc_ptrs.size();
    for (int c = 0; c < n_vcs; ++c)
    if (this->vc_ptrs[c] != nullptr)
    {
        delete this->vc_ptrs[c];
    }
    this->vc_ptrs.clear();
    this->ids.clear();
    this->seeds.clear();
}
// ====================================================================

// AUXILIARY METHODS ==================================================
int Mesh::get_new_id() const
{
    // PARAMETERS -------------------------
    const int n_cells = this->cells.size();
    // ------------------------------------

    // VARIABLES
    int id = 0;
    // ---------

    // LOOP OVER THE CELLS ------------------
    for (int c = 0; c < n_cells; ++c)
    {
        id = std::max(id, this->cells[c].id);
    }
    // --------------------------------------

    return (id+1);
}

int Mesh::find_cell_index_by_id(const int id) const
{
    // PARAMETERS -------------------------
    const int n_cells = this->cells.size();
    // ------------------------------------

    // SEARCH FOR THE CELL WITH MATCHING ID
    int c = 0;
    bool found = false;

    while (!found && c < n_cells)
    {
        if (this->cells[c].id == id)
        {
            found = true;
        }
        else
        {
            c += 1;
        }
    }
    // ------------------------------------

    return c;
}

int Mesh::find_wall_index_by_id(const int id) const
{
    // PARAMETERS -------------------------
    const int n_walls = this->walls.size();
    // ------------------------------------

    // SEARCH FOR THE WALL WITH MATCHING ID
    int w = 0;
    bool found = false;

    while (!found && w < n_walls)
    {
        if (this->walls[w].id == id)
        {
            found = true;
        }
        else
        {
            w += 1;
        }
    }
    // ------------------------------------

    return w;
}
// ====================================================================

// INITIALIZATION: TESSELLATION =======================================
void Mesh::init_tessellation()
{
    // FREE PREVIOUS ONE ------
    this->clear_tessellation();
    // ------------------------

    // PARAMETERS -------------------------
    const int n_vcs = this->vc_ptrs.size();
    // ------------------------------------

    // VORONOI CELLS --------------------------------------------------
    for (int k = 0; k < n_vcs; ++k)
    {
        voro::voronoicell_neighbor * vc = this->vc_ptrs[k];

        // CREATE A VORONOI CELL
        VoronoiCell cell;

        // Set the id
        cell.id = this->ids[k];

        if (vc != nullptr)
        {
            // Set the neighbors
            vc->neighbors(cell.nbr);

            // Compute the centroid
            vc->centroid(cell.centroid[0], cell.centroid[1], cell.centroid[2]);
            cell.centroid[0] += this->seeds[3*k+0];
            cell.centroid[1] += this->seeds[3*k+1];
            cell.centroid[2] += this->seeds[3*k+2];

            // Compute the volume
            cell.volume = vc->volume();
        }

        // ADD THE VORONOI CELL
        this->cells.push_back(cell);
    }
    // ----------------------------------------------------------------

    // FIX NEIGHBORS CONNECTIVITY -------------------------------------
    for (int c = 0; c < n_vcs; ++c)
    {
        VoronoiCell & cell = this->cells[c];
        const int n_nbr = cell.nbr.size();

        for (int n = 0; n < n_nbr; ++n)
        {
            int nbr_id = cell.nbr[n];

            // WALL
            if (nbr_id < 0)
            {

            }
            // THE CELL IS NEIGHBOR WITH ITSELF
            else if (nbr_id == cell.id)
            {
                io::error("mesh.cpp - Mesh::init_tessellation",
                          "Cells neighboring with themselves must be handled yet.");
            }
            // AN ACTUAL NEIGHBORING CELL
            else
            {
                int nbr_c = this->find_cell_index_by_id(nbr_id);

                if (nbr_c == n_vcs)
                {
                    std::string msg = "The neighbor cells information is not consistent:\n";
                    msg += "| current cell id: "+std::to_string(cell.id)+";\n";
                    msg += "| neighbor cell id: "+std::to_string(nbr_id)+".\n";
                    io::error("mesh.cpp - Mesh::init_tessellation", msg);
                }

                cell.nbr[n] = nbr_c;
            }
        }
    }
    // ----------------------------------------------------------------
    
    // VORONOI FACES --------------------------------------------------
    for (int c = 0; c < n_vcs; ++c)
    {
        const double * S = &this->seeds[3*c];
        voro::voronoicell_neighbor * vc = this->vc_ptrs[c];
        VoronoiCell & cell = this->cells[c];
        const int n_nbr = cell.nbr.size();
        const int n_cell_faces = n_nbr;

        if (vc != nullptr)
        {
            std::vector<double> vertices;
            std::vector<int> face_vertices;

            vc->vertices(S[0], S[1], S[2], vertices);
            vc->face_vertices(face_vertices);

            cell.f_conn.resize(n_cell_faces);

            int face_vertices_pos = 0;
            for (int cf = 0; cf < n_cell_faces; ++cf)
            {
                const int nbr_c = cell.nbr[cf];
                const bool face_is_wall = (nbr_c < 0);
                const bool new_face = (face_is_wall || (nbr_c > c));
                const int n_face_vertices = face_vertices[face_vertices_pos];
                const int bou_type = ((face_is_wall) ? BOU_TYPE_WALL : BOU_TYPE_INTERFACE);

                if (new_face)
                {
                    // CREATE A FACE
                    VoronoiFace face;

                    // Set the boundary type
                    face.bou_type = bou_type;

                    // Set the parent cells
                    face.parent_cells[0] = c;
                    face.parent_cells[1] = nbr_c;

                    // Eval the centroid and the area
                    face.centroid[0] = 0.0;
                    face.centroid[1] = 0.0;
                    face.centroid[2] = 0.0;
                    face.area = 0.0;

                    const int v0 = face_vertices[face_vertices_pos+1];
                    const double * V0 = &vertices[3*v0];

                    for (int fv = 1; fv < (n_face_vertices-1); ++fv)
                    {
                        const int vi = face_vertices[face_vertices_pos+1+fv];
                        const int vj = face_vertices[face_vertices_pos+1+fv+1];
                        const double * Vi = &vertices[3*vi];
                        const double * Vj = &vertices[3*vj];
                        const double tri_area = math::tri_area_3d(V0, Vi, Vj);
                        double tri_C[3];
                        math::tri_centroid(V0, Vi, Vj, tri_C);

                        face.centroid[0] += tri_C[0]*tri_area;
                        face.centroid[1] += tri_C[1]*tri_area;
                        face.centroid[2] += tri_C[2]*tri_area;

                        face.area += tri_area;
                    }

                    face.centroid[0] *= 1.0/face.area;
                    face.centroid[1] *= 1.0/face.area;
                    face.centroid[2] *= 1.0/face.area;

                    // Store the index in the cell's face connectivity
                    cell.f_conn[cf] = this->faces.size();

                    // ADD THE VORONOI FACE
                    this->faces.push_back(face);
                }
                else
                {
                    // FIND THE FACE
                    const int n_faces = this->faces.size();

                    bool found = false;
                    int ff = 0;
                    while ((!found) && (ff < n_faces))
                    {
                        const VoronoiFace & face = this->faces[ff];

                        if ((face.parent_cells[0] == nbr_c) &&
                            (face.parent_cells[1] == c) &&
                            (face.bou_type != BOU_TYPE_WALL))
                        {
                            found = true;

                            // Store the index in the cell's face
                            // connectivity
                            cell.f_conn[cf] = ff;
                        }
                        else
                        {
                            ff += 1;
                        }
                    }

                    if (!found)
                    {
                        std::string msg = "The neighbor cells information is not consistent:\n";
                        msg += "| current cell:\n";
                        msg += "| - id(index): "+std::to_string(cell.id)+"("+std::to_string(c)+");\n";
                        msg += "| - nbr:"; for (int i = 0; i < n_nbr; ++i) msg += " "+std::to_string(cell.nbr[i])+((i == cf) ? "*" : ""); msg += ";\n";
                        msg += "| - faces:"; for (int i = 0; i < n_cell_faces; ++i) msg += " "+std::to_string(cell.f_conn[i])+((i == cf) ? "*" : ""); msg += ";\n";
                        
                        const VoronoiCell & nbr_cell = this->cells[nbr_c];
                        const int nbr_n_nbr = nbr_cell.nbr.size();
                        const int nbr_n_cell_faces = nbr_cell.f_conn.size();
                        msg += "| neighbor cell:\n";
                        msg += "| - id(index): "+std::to_string(nbr_cell.id)+"("+std::to_string(nbr_c)+");\n";
                        msg += "| - nbr:"; for (int i = 0; i < nbr_n_nbr; ++i) msg += " "+std::to_string(nbr_cell.nbr[i]); msg += ";\n";
                        msg += "| - faces:"; for (int i = 0; i < nbr_n_cell_faces; ++i) msg += " "+std::to_string(nbr_cell.f_conn[i]); msg += ";\n";

                        msg += "| current face:\n";
                        msg += "| - index: "+std::to_string(cf)+";\n";
                        io::error("mesh.cpp - Mesh::init_tessellation",
                                   msg);
                    }
                }

                // MOVE TO THE NEXT FACE
                face_vertices_pos += n_face_vertices+1;
            }
        }
    }
    // ----------------------------------------------------------------
    
    // VORONOI WALLS --------------------------------------------------
    {
        const int n_faces = this->faces.size();

        for (int f = 0; f < n_faces; ++f)
        {
            const VoronoiFace & face = this->faces[f];

            if (face.bou_type == BOU_TYPE_WALL)
            {
                const int wall_id = -face.parent_cells[1];

                if (wall_id < 0)
                {
                    io::error("mesh.cpp - Mesh::init_tessellation",
                              "A wall face must have a negative wall id associated to it.");
                }

                // CHECK WHETHER WE HAVE ADDED A WALL WITH THIS ID
                const size_t w = this->find_wall_index_by_id(wall_id);

                if (w == this->walls.size())
                {
                    // CREATE A WALL
                    VoronoiWall wall;

                    // Set the id
                    wall.id = wall_id;

                    // Add the face
                    wall.f_conn.clear();
                    wall.f_conn.push_back(f);
                    
                    // Add the cell
                    wall.c_conn.clear();
                    wall.c_conn.push_back(face.parent_cells[0]);

                    // Push the wall
                    this->walls.push_back(wall);
                }
                else
                {
                    // GET THE WALL
                    VoronoiWall & wall = this->walls[w];

                    // Add the face
                    wall.f_conn.push_back(f);
                    
                    // Add the cell
                    wall.c_conn.push_back(face.parent_cells[0]);
                }
            }
        }
    }
    // ----------------------------------------------------------------
}
// ====================================================================

// EXPORT TO VTK FORMAT: TESSELLATION =================================
void Mesh::export_tessellation_vtk(const std::string & filepath) const
{
    // PARAMETERS -----------------------------------------------------
    const int n_vcs = this->vc_ptrs.size();

    const std::vector<std::string> VTK_cell_fields_name = 
    {
        "DIM",
        "BOU",
        "PARENT",
        "PARENT_GROUP"
    };
    const int n_cell_fields = VTK_cell_fields_name.size();
    // ----------------------------------------------------------------

    // VARIABLES ------------------------------------------------------
    vtk::Cell_conn_t n_VTK_points;
    vtk::Cell_offs_t n_VTK_cells;

    long VTK_cell_conn_len;

    vtk::Cell_offs_t VTK_cells_pos;

    std::vector<vtk::Float_t> VTK_points;
    std::vector<vtk::Cell_conn_t> VTK_cell_conn;
    std::vector<vtk::Cell_offs_t> VTK_cell_offset;
    std::vector<vtk::Cell_type_t> VTK_cell_type;
    std::vector<std::vector<vtk::Int_t>> VTK_cell_fields(n_cell_fields);

    int af;
    // ----------------------------------------------------------------

    // VTK POINTS AND VTK CELLS ---------------------------------------
    n_VTK_points = 0;
    n_VTK_cells = 0;

    VTK_cell_offset.push_back(0);

    af = 0;

    for (int c = 0; c < n_vcs; ++c)
    {
        const int id = this->ids[c];

        if (this->vc_ptrs[c] != nullptr)
        {
            const double * S = &this->seeds[3*c];
            voro::voronoicell_neighbor * vc = this->vc_ptrs[c];
            
            std::vector<double> vertices;
            std::vector<int> face_vertices;
            std::vector<int> nbr;
            int n_vertices, n_faces;

            vc->vertices(S[0], S[1], S[2], vertices);
            vc->face_vertices(face_vertices);
            vc->neighbors(nbr);
            n_vertices = vertices.size()/3;
            n_faces = nbr.size();

            // VTK CELLS: ADD FACE TRIANGLES
            int face_vertices_pos = 0;
            for (int f = 0; f < n_faces; ++f)
            {
                const int nbr_id = nbr[f];
                const bool face_is_wall = (nbr_id < 0);
                const bool new_face = (face_is_wall || (nbr_id > id));
                const int n_face_vertices = face_vertices[face_vertices_pos];
                const int bou_type = this->faces[this->cells[c].f_conn[f]].bou_type;

                if (new_face)
                {
                    const int v0 = face_vertices[face_vertices_pos+1];

                    for (int fv = 1; fv < (n_face_vertices-1); ++fv)
                    {
                        const int vi = face_vertices[face_vertices_pos+1+fv];
                        const int vj = face_vertices[face_vertices_pos+1+fv+1];

                        VTK_cell_conn.push_back(n_VTK_points+v0);
                        VTK_cell_conn.push_back(n_VTK_points+vi);
                        VTK_cell_conn.push_back(n_VTK_points+vj);

                        VTK_cell_offset.push_back(VTK_cell_offset.back()+3);

                        VTK_cell_type.push_back(VTK_TRIANGLE);

                        VTK_cell_fields[0].push_back(2);
                        VTK_cell_fields[1].push_back(bou_type);
                        VTK_cell_fields[2].push_back(af);

                        if (bou_type == BOU_TYPE_WALL)
                        {
                            VTK_cell_fields[3].push_back(-nbr_id);
                        }
                        else if (bou_type == BOU_TYPE_INTERFACE)
                        {
                            VTK_cell_fields[3].push_back(0);
                        }
                        else if (bou_type == BOU_TYPE_CRACK)
                        {
                            const int n_cracks = this->cracks.size();
                            std::vector<int>::const_iterator it;
                            int cr = 0;
                            bool found = false;
                            while ((!found) && (cr < n_cracks))
                            {
                                it = std::find(this->cracks[cr].f_conn.begin(),
                                               this->cracks[cr].f_conn.end(),
                                               this->cells[c].f_conn[f]);
                                if (it != this->cracks[cr].f_conn.end())
                                {
                                    found = true;
                                }
                                else
                                {
                                    cr += 1;
                                }
                            }

                            if (cr < n_cracks)
                            {
                                VTK_cell_fields[3].push_back(this->cracks[cr].id);
                            }
                            else
                            {
                                VTK_cell_fields[3].push_back(-1);
                            }
                        }
                        else
                        {
                            VTK_cell_fields[3].push_back(-1);
                        }

                        n_VTK_cells += 1;
                    }

                    af += 1;
                }

                // MOVE TO THE NEXT FACE
                face_vertices_pos += n_face_vertices+1;
            }

            // VTK CELLS: ADD CELL TETRAHEDRA
            {
                int cnt;
                const int v0 = 0;

                face_vertices_pos = 0;
                for (int f = 0; f < n_faces; ++f)
                {
                    const int n_face_vertices = face_vertices[face_vertices_pos];

                    cnt = std::count(&face_vertices[face_vertices_pos+1],
                                     &face_vertices[face_vertices_pos+n_face_vertices+1],
                                     v0);

                    if (cnt == 0)
                    {
                        const int vi = face_vertices[face_vertices_pos+1];

                        for (int fv = 1; fv < (n_face_vertices-1); ++fv)
                        {
                            const int vj = face_vertices[face_vertices_pos+1+fv];
                            const int vk = face_vertices[face_vertices_pos+1+fv+1];

                            {
                                VTK_cell_conn.push_back(n_VTK_points+v0);
                                VTK_cell_conn.push_back(n_VTK_points+vi);
                                VTK_cell_conn.push_back(n_VTK_points+vj);
                                VTK_cell_conn.push_back(n_VTK_points+vk);

                                VTK_cell_offset.push_back(VTK_cell_offset.back()+4);

                                VTK_cell_type.push_back(VTK_TETRA);

                                VTK_cell_fields[0].push_back(3);
                                VTK_cell_fields[1].push_back(-1);
                                VTK_cell_fields[2].push_back(id);
                                VTK_cell_fields[3].push_back(0);

                                n_VTK_cells += 1;
                            }
                        }

                        af += 1;
                    }

                    // MOVE TO THE NEXT FACE
                    face_vertices_pos += n_face_vertices+1;
                }
            }

            // VTK POINT
            n_VTK_points += n_vertices;

            for (int v = 0; v < n_vertices; ++v)
            {
                VTK_points.push_back(vertices[3*v+0]);
                VTK_points.push_back(vertices[3*v+1]);
                VTK_points.push_back(vertices[3*v+2]);
            }
        }
    }
    // ----------------------------------------------------------------

    // DUMP DATA ------------------------------------------------------
    vtk::print_unstructured_data(filepath,
                                 n_VTK_points,
                                 n_VTK_cells,
                                 VTK_points,
                                 VTK_cell_conn,
                                 VTK_cell_offset,
                                 VTK_cell_type,
                                 VTK_cell_fields,
                                 VTK_cell_fields_name,
                                 {},
                                 {},
                                 "ascii");
    // ----------------------------------------------------------------
}
// ====================================================================
// ####################################################################

}