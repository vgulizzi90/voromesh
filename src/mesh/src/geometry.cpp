// geometry.cpp

#include <limits>
#include "geometry.hpp"
#include "utils_vtk.hpp"

namespace voromesh
{
// VORONOI CELL FACE CLASS ############################################
// DESTRUCTOR =========================================================
VorocellFace::~VorocellFace()
{

}
// ====================================================================
// ####################################################################



// VORONOI CELL CLASS #################################################
// DESTRUCTOR =========================================================
Vorocell::~Vorocell()
{
    if (this->pp != nullptr) delete this->pp;
}
// ====================================================================

// EVAL PROPERTIES USING VORO++ METHODS ===============================
/*
void Vorocell::eval_properties()
{
    if (this->pp != nullptr)
    {
        this->pp->vertices(this->centroid[0], this->centroid[1], this->centroid[2],
                           this->vertices);

        this->pp->face_vertices(this->face_vertices);

        this->pp->neighbors(this->nbr);
    }
}
*/
// ====================================================================

// EVAL SHORTEST EDGE =================================================
/*
double Vorocell::eval_shortest_edge()
{
    // VARIABLES
    const double * A, * B;
    double l;
    double shortest_edge = std::numeric_limits<double>::max();

    if (this->pp != nullptr)
    {
        for (int i = 0; i < this->pp->p; ++i)
        {
            for (int j = 0; j < this->pp->nu[i]; ++j)
            {
                if ((this->pp->ed[i][j]) > i)
                {

                }
                else
                {
                    A = &this->vertices[3*i];
                    B = &this->vertices[3*(this->pp->ed[i][j])];
                    l = math::L2_distance(A, B, 3);
                    shortest_edge = std::min(shortest_edge, l);
                }
            }
        }
    }

    return shortest_edge;
}
*/
// ====================================================================
// ####################################################################


// GEOMETRY CLASS #####################################################
// CONSTRUCTOR ========================================================
Geometry::Geometry()
{

}
// ====================================================================

// DESTRUCTOR =========================================================
Geometry::~Geometry()
{
    this->free();
}
// ====================================================================

// FREE MEMORY ========================================================
void Geometry::free()
{
    // CLEAR ALL VORONOI CELLS AND CELL FACES
    this->faces.clear();
    this->cells.clear();
}
// ====================================================================

// INITIALIZATION =====================================================
// ====================================================================

// BUILD THE GEOMETRICAL ENTITIES =====================================
void Geometry::build()
{
    // PARAMETER
    const int n_vc = this->cells.size();

    // VARIABLES
    
    // FREE MEMORY
    this->faces.clear();

    // EVAL THE CENTROIDS OF THE VORONOI CELLS ------------------------
    for (int k = 0; k < n_vc; ++k)
    {
        Vorocell & vc = this->cells[k];

        if (vc.pp != nullptr)
        {
            vc.pp->centroid(vc.centroid[0], vc.centroid[1], vc.centroid[2]);
            vc.centroid[0] += vc.seed[0];
            vc.centroid[1] += vc.seed[1];
            vc.centroid[2] += vc.seed[2];
        }
    }
    // ----------------------------------------------------------------

    // EVAL THE NEIGHBORS ---------------------------------------------
    for (int k = 0; k < n_vc; ++k)
    {
        Vorocell & vc = this->cells[k];
        
        if (vc.pp != nullptr)
        {
            vc.pp->neighbors(vc.nbr);

            // MAKE SURE THE CONNECTIVITY IS CONSISTENT
            const int n_nbr = vc.nbr.size();

            for (int n = 0; n < n_nbr; ++n)
            {
                int nbr_id = vc.nbr[n];
                // WALL
                if (nbr_id < 0)
                {

                }
                // THE CELL IS NEIGHBOR WITH ITSELF
                else if (nbr_id == vc.id)
                {
                    io::error("geometry.cpp - Geometry::build",
                              "Cells neighboring with themselves must be handled yet.");
                }
                // AN ACTUAL NEIGHBORING CELL
                else
                {
                    int nbr_k = 0;
                    bool found = false;
                    while ((!found) && (nbr_k < n_vc))
                    {
                        Vorocell & nbr_vc = this->cells[nbr_k];
                        if (nbr_vc.id == nbr_id)
                        {
                            found = true;
                            vc.nbr[n] = nbr_k;
                        }
                        else
                        {
                            nbr_k += 1;
                        }
                    }

                    if (!found)
                    {
                        io::error("geometry.cpp - Geometry::build",
                                  "The neighbor cells information is not consistent.");
                    }
                }
            }
        }
    }
    // ----------------------------------------------------------------

    // ADD THE VORONOI CELL FACES -------------------------------------
    for (int k = 0; k < n_vc; ++k)
    {
        Vorocell & vc = this->cells[k];
        
        if (vc.pp != nullptr)
        {
            const int n_cell_faces = vc.nbr.size();

            std::vector<double> ver;
            std::vector<int> face_vertices;
            int face_vertice_pos;
            vc.pp->vertices(vc.seed[0], vc.seed[1], vc.seed[2], ver);
            vc.pp->face_vertices(face_vertices);

            vc.f_conn.resize(n_cell_faces);
            vc.f_ori.resize(n_cell_faces);

            face_vertice_pos = 0;
            for (int f = 0; f < n_cell_faces; ++f)
            {
                const int nbr_k = vc.nbr[f];
                const bool face_is_wall = (nbr_k < 0);
                const bool new_face = (face_is_wall || (nbr_k > k));

                // ADD A NEW FACE
                if (new_face)
                {
                    // CREATE A FACE
                    VorocellFace vcf;

                    // FACE INDEX
                    vc.f_conn[f] = this->faces.size();

                    // FACE VERTICES
                    const int n_face_vertices = face_vertices[face_vertice_pos];
                    vcf.vertices.resize(3*n_face_vertices);
                    for (int v = 0; v < n_face_vertices; ++v)
                    {
                        const int vv = face_vertices[face_vertice_pos+1+v];
                        vcf.vertices[3*v+0] = ver[3*vv+0];
                        vcf.vertices[3*v+1] = ver[3*vv+1];
                        vcf.vertices[3*v+2] = ver[3*vv+2];
                    }
                    
                    // FIRST AND SECOND VERTICES OF THE FACE
                    const double * V0 = &vcf.vertices[3*0];
                    const double * V1 = &vcf.vertices[3*1];

                    // COMPUTE FACE CENTROID
                    double TC[3];
                    double face_area = 0.0;

                    vcf.centroid[0] = 0.0;
                    vcf.centroid[1] = 0.0;
                    vcf.centroid[2] = 0.0;

                    for (int t = 0; t < (n_face_vertices-2); ++t)
                    {
                        double tri_area;
                        const double * Vi = &vcf.vertices[3*(t+1)];
                        const double * Vj = &vcf.vertices[3*(t+2)];

                        math::tri_centroid(V0, Vi, Vj, TC);
                        tri_area = math::tri_area_3d(V0, Vi, Vj);

                        face_area += tri_area;
                        vcf.centroid[0] += TC[0]*tri_area;
                        vcf.centroid[1] += TC[1]*tri_area;
                        vcf.centroid[2] += TC[2]*tri_area;
                    }
                    vcf.centroid[0] /= face_area;
                    vcf.centroid[1] /= face_area;
                    vcf.centroid[2] /= face_area;

                    // PARENT CELLS INFO
                    vcf.parent_cells[0] = k;
                    vcf.parent_cells[1] = nbr_k;

                    // BOUNDARY TYPE
                    vcf.bou_type = ((face_is_wall)? BOU_TYPE_WALL : BOU_TYPE_INTERFACE);

                    // FACE ORIENTATION
                    double tmp, tmp_V[3];
                    math::cross(V0, V1, vcf.centroid, tmp_V);
                    tmp  = tmp_V[0]*(vcf.centroid[0]-vc.centroid[0]);
                    tmp += tmp_V[1]*(vcf.centroid[1]-vc.centroid[1]);
                    tmp += tmp_V[2]*(vcf.centroid[2]-vc.centroid[2]);

                    const bool change_vertices_order = (tmp < 0.0);

                    if (change_vertices_order)
                    {
                        std::vector<double> aux_face_ver(vcf.vertices);
                        for (int v = 0; v < n_face_vertices; ++v)
                        {
                            const int vv = n_face_vertices-1-v;
                            vcf.vertices[3*v+0] = aux_face_ver[3*vv+0];
                            vcf.vertices[3*v+1] = aux_face_ver[3*vv+1];
                            vcf.vertices[3*v+2] = aux_face_ver[3*vv+2];
                        }
                    }
                    
                    // ADD THE FACE
                    this->faces.push_back(vcf);
                }
                // THE FACE ALREADY EXISTS
                else
                {
                    // FIND THE FACE
                    const int n_faces = this->faces.size();
                    bool found = false;
                    int ff = 0;
                    while ((!found) && (ff < n_faces))
                    {
                        const VorocellFace & vcf = this->faces[ff];

                        if ((vcf.parent_cells[0] == nbr_k) &&
                            (vcf.parent_cells[1] == k) &&
                            (vcf.bou_type != 0))
                        {
                            found = true;

                            // CELL INFORMATION
                            vc.f_conn[f] = ff;
                        }
                        else
                        {
                            ff += 1;
                        }
                    }

                    if (!found)
                    {
                        io::error("geometry.cpp - Geometry::build",
                                  "The neighbor cells information is not consistent.");
                    }
                }

                // MOVE TO THE NEXT FACE
                face_vertice_pos += face_vertices[face_vertice_pos]+1;
            }
        }
    }
    // ----------------------------------------------------------------
}
// ====================================================================

// EVAL THE GEOMETRIC TOLERANCE =======================================
double Geometry::eval_tolerance(const double rel_tol)
{
    // PARAMETERS
    const int n_vc = this->cells.size();
    
    // VARIABLES
    double shortest_edge = std::numeric_limits<double>::max();

    // [...]
    io::error("geometry.cpp", "to be implemented");

    // RETURN
    return shortest_edge*rel_tol;
}
// ====================================================================

// EXPORT TO VTK FORMAT ===============================================
void Geometry::export_VTK(const std::string & filepath) const
{
    // PARAMETERS -----------------------------------------------------
    const int n_vc = this->cells.size();
    const int n_vcf = this->faces.size();

    const std::vector<std::string> VTK_cell_field_names = {"ID", "DIM_TYPE", "N_NBR"};
    const int n_cell_fields = VTK_cell_field_names.size();
    // ----------------------------------------------------------------

    // VARIABLES ------------------------------------------------------
    vtk::Cell_conn_t VTK_nodes_pos;
    vtk::Cell_offs_t VTK_cells_pos;
    vtk::Cell_conn_t n_VTK_nodes;
    vtk::Cell_offs_t n_VTK_cells;
    long VTK_cell_conn_len;

    std::vector<vtk::Float_t> VTK_nodes;
    std::vector<vtk::Cell_conn_t> VTK_cell_conn;
    std::vector<vtk::Cell_offs_t> VTK_cell_offset;
    std::vector<vtk::Cell_type_t> VTK_cell_type;
    std::vector<std::vector<vtk::Int_t>> VTK_cell_field(n_cell_fields);

    std::vector<long> aux_face_conn(this->faces.size(), 0);
    // ----------------------------------------------------------------

    // COUNT THE NUMBER OF VTK NODES ----------------------------------
    n_VTK_nodes = 0;

    for (int f = 0; f < n_vcf; ++f)
    {
        const VorocellFace & vcf = this->faces[f];
        n_VTK_nodes += vcf.vertices.size()/3;
    }

    n_VTK_nodes += n_vc;
    // ----------------------------------------------------------------

    // COUNT THE NUMBER OF VTK NODES AND CELLS ------------------------
    n_VTK_cells = 0;
    VTK_cell_conn_len = 0L;

    for (int f = 0; f < n_vcf; ++f)
    {
        const VorocellFace & vcf = this->faces[f];
        const int n_face_vertices = vcf.vertices.size()/3;
        n_VTK_cells += n_face_vertices-2;
        VTK_cell_conn_len += 3*(n_face_vertices-2);
    }

    for (int c = 0; c < n_vc; ++c)
    {
        const Vorocell & vc = this->cells[c];
        const int n_cell_faces = vc.f_conn.size();
        for (int cf = 0; cf < n_cell_faces; ++cf)
        {
            const int f = vc.f_conn[cf];
            const VorocellFace & vcf = this->faces[f];
            const int n_face_vertices = vcf.vertices.size()/3;

            n_VTK_cells += n_face_vertices-2;
            VTK_cell_conn_len += 4*(n_face_vertices-2);
        }
    }
    // ----------------------------------------------------------------

    // RESIZE MEMORY --------------------------------------------------
    VTK_nodes.resize(3*n_VTK_nodes);
    VTK_cell_conn.resize(VTK_cell_conn_len);
    VTK_cell_offset.resize(n_VTK_cells+1);
    VTK_cell_type.resize(n_VTK_cells);
    for (int f = 0; f < n_cell_fields; ++f)
    {
        VTK_cell_field[f].resize(n_VTK_cells);
    }
    // ----------------------------------------------------------------

    // EVAL NODES -----------------------------------------------------
    vtk::Float_t * VTK_nodes_ptr = VTK_nodes.data();

    for (int f = 0; f < n_vcf; ++f)
    {
        const VorocellFace & vcf = this->faces[f];
        const int n_face_vertices = vcf.vertices.size()/3;

        for (int fv = 0; fv < n_face_vertices; ++fv)
        {
            *VTK_nodes_ptr = vcf.vertices[3*fv+0]; ++VTK_nodes_ptr;
            *VTK_nodes_ptr = vcf.vertices[3*fv+1]; ++VTK_nodes_ptr;
            *VTK_nodes_ptr = vcf.vertices[3*fv+2]; ++VTK_nodes_ptr;
        }
    }

    for (int c = 0; c < n_vc; ++c)
    {
        const Vorocell & vc = this->cells[c];

        *VTK_nodes_ptr = vc.centroid[0]; ++VTK_nodes_ptr;
        *VTK_nodes_ptr = vc.centroid[1]; ++VTK_nodes_ptr;
        *VTK_nodes_ptr = vc.centroid[2]; ++VTK_nodes_ptr;
    }
    // ----------------------------------------------------------------

    // EVAL CELLS -----------------------------------------------------
    vtk::Cell_conn_t * VTK_cell_conn_ptr = VTK_cell_conn.data();
    vtk::Cell_offs_t * VTK_cell_offset_ptr = VTK_cell_offset.data();
    *VTK_cell_offset_ptr = 0; ++VTK_cell_offset_ptr;
    vtk::Cell_type_t * VTK_cell_type_ptr = VTK_cell_type.data();
    
    std::vector<vtk::Int_t *> VTK_cell_field_ptr(n_cell_fields);
    for (int fld = 0; fld < n_cell_fields; ++fld)
    {
        VTK_cell_field_ptr[fld] = VTK_cell_field[fld].data();
    }

    VTK_nodes_pos = 0;
    VTK_cells_pos = 0;
    for (int f = 0; f < n_vcf; ++f)
    {
        const VorocellFace & vcf = this->faces[f];
        const int n_face_vertices = vcf.vertices.size()/3;

        for (int fv = 0; fv < (n_face_vertices-2); ++fv)
        {
            *VTK_cell_conn_ptr = VTK_nodes_pos+0; ++VTK_cell_conn_ptr;
            *VTK_cell_conn_ptr = VTK_nodes_pos+fv+1; ++VTK_cell_conn_ptr;
            *VTK_cell_conn_ptr = VTK_nodes_pos+fv+2; ++VTK_cell_conn_ptr;

            VTK_cells_pos += 3;
            *VTK_cell_offset_ptr = VTK_cells_pos; ++VTK_cell_offset_ptr;

            *VTK_cell_type_ptr = 5; ++VTK_cell_type_ptr;

            *(VTK_cell_field_ptr[0]) = f; ++(VTK_cell_field_ptr[0]);
            *(VTK_cell_field_ptr[1]) = 2; ++(VTK_cell_field_ptr[1]);
            *(VTK_cell_field_ptr[2]) = 2; ++(VTK_cell_field_ptr[2]);
        }

        aux_face_conn[f] = VTK_nodes_pos;

        VTK_nodes_pos += n_face_vertices;
    }
    
    for (int c = 0; c < n_vc; ++c)
    {
        const Vorocell & vc = this->cells[c];
        const int n_cell_faces = vc.f_conn.size();
        for (int cf = 0; cf < n_cell_faces; ++cf)
        {
            const int f = vc.f_conn[cf];
            const VorocellFace & vcf = this->faces[f];
            const int n_face_vertices = vcf.vertices.size()/3;

            for (int fv = 0; fv < (n_face_vertices-2); ++fv)
            {
                *VTK_cell_conn_ptr = VTK_nodes_pos; ++VTK_cell_conn_ptr;
                *VTK_cell_conn_ptr = aux_face_conn[f]+0; ++VTK_cell_conn_ptr;
                *VTK_cell_conn_ptr = aux_face_conn[f]+fv+1; ++VTK_cell_conn_ptr;
                *VTK_cell_conn_ptr = aux_face_conn[f]+fv+2; ++VTK_cell_conn_ptr;

                VTK_cells_pos += 4;
                *VTK_cell_offset_ptr = VTK_cells_pos; ++VTK_cell_offset_ptr;

                *VTK_cell_type_ptr = 10; ++VTK_cell_type_ptr;

                *(VTK_cell_field_ptr[0]) = vc.id; ++(VTK_cell_field_ptr[0]);
                *(VTK_cell_field_ptr[1]) = 3; ++(VTK_cell_field_ptr[1]);
                *(VTK_cell_field_ptr[2]) = n_cell_faces; ++(VTK_cell_field_ptr[2]);
            }
        }

        VTK_nodes_pos += 1;
    }
    // ----------------------------------------------------------------

    // DUMP DATA ------------------------------------------------------
    vtk::print_unstructured_data(filepath,
                                 n_VTK_nodes,
                                 n_VTK_cells,
                                 VTK_nodes,
                                 VTK_cell_conn,
                                 VTK_cell_offset,
                                 VTK_cell_type,
                                 VTK_cell_field,
                                 VTK_cell_field_names,
                                 {},
                                 {},
                                 "ascii");
    // ----------------------------------------------------------------
}
// ====================================================================
// ####################################################################


// GEOMETRY ROUTINES ##################################################
// ####################################################################

}