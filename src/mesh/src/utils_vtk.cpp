// utils_vtk.cpp

#include "utils_vtk.hpp"

namespace voromesh
{
namespace vtk
{

void print_unstructured_data(const std::string & filepath,
                             const Cell_conn_t & n_nodes,
                             const Cell_offs_t & n_cells,
                             const std::vector<Float_t> & nodes,
                             const std::vector<Cell_conn_t> & cell_conn,
                             const std::vector<Cell_offs_t> & cell_offset,
                             const std::vector<Cell_type_t> & cell_type,
                             const std::vector<std::vector<Int_t>> & cell_field,
                             const std::vector<std::string> & cell_field_name,
                             const std::vector<std::vector<Float_t>> & nodal_field,
                             const std::vector<std::string> & nodal_field_name,
                             const std::string & fmt)
{
    // VARIABLES ------------------------------------------------------
    std::ofstream fp;
    // ----------------------------------------------------------------

    // CHECK FORMAT ---------------------------------------------------
    if ((fmt.compare("ascii") != 0) && (fmt.compare("binary") != 0))
    {
        io::error("utils_vtk.cpp - print_unstructured_data",
                  "Unexpected file format: "+fmt);
    }
    // ----------------------------------------------------------------

    // OPEN FILE ------------------------------------------------------
    fp.open(filepath.c_str(), std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);

    if (!fp.good())
    {
        io::error("utils_vtk.cpp - print_unstructured_data",
                  "Could not open file: "+filepath);
    }
    // ----------------------------------------------------------------

    // ASCII FORMAT ---------------------------------------------------
    if (fmt.compare("ascii") == 0)
    {
        print_unstructured_data_ascii(fp,
                                      n_nodes,
                                      n_cells,
                                      nodes,
                                      cell_conn,
                                      cell_offset,
                                      cell_type,
                                      cell_field,
                                      cell_field_name,
                                      nodal_field,
                                      nodal_field_name);
    }
    // ----------------------------------------------------------------
    // BINARY FORMAT --------------------------------------------------
    else if (fmt.compare("binary") == 0)
    {
        print_unstructured_data_binary(fp,
                                       n_nodes,
                                       n_cells,
                                       nodes,
                                       cell_conn,
                                       cell_offset,
                                       cell_type,
                                       cell_field,
                                       cell_field_name,
                                       nodal_field,
                                       nodal_field_name);
    }
    // ----------------------------------------------------------------

    // CLOSE FILE
    fp.close();
    // ----------
}

void print_unstructured_data_ascii(std::ofstream & fp,
                                   const Cell_conn_t & n_nodes,
                                   const Cell_offs_t & n_cells,
                                   const std::vector<Float_t> & nodes,
                                   const std::vector<Cell_conn_t> & cell_conn,
                                   const std::vector<Cell_offs_t> & cell_offset,
                                   const std::vector<Cell_type_t> & cell_type,
                                   const std::vector<std::vector<Int_t>> & cell_field,
                                   const std::vector<std::string> & cell_field_name,
                                   const std::vector<std::vector<Float_t>> & nodal_field,
                                   const std::vector<std::string> & nodal_field_name)
{
    // PARAMETERS -----------------------------------------------------
    const int cell_offset_stride = 20;
    const int cell_type_stride = 20;
    const int cell_field_stride = 20;
    const int nodal_field_stride = 10;

    const int n_cell_fields = cell_field.size();
    const int n_nodal_fields = nodal_field.size();
    // ----------------------------------------------------------------

    // SETTINGS -------------------------------------------------------
    fp.precision(17);
    // ----------------------------------------------------------------

    // HEADER ---------------------------------------------------------
    fp << "<?xml version=\"1.0\"?>\n";
    fp << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    fp << "<UnstructuredGrid>\n";

    fp << "<Piece NumberOfPoints=\"" << n_nodes << "\" NumberOfCells=\"" << n_cells << "\">\n";
    // ----------------------------------------------------------------

    // NODES ----------------------------------------------------------
    fp << "<Points>\n";
    fp << "  <DataArray type=\"" << Float_t_description << "\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (Cell_conn_t n = 0; n < n_nodes; ++n)
    {
        fp << "  " << nodes[3*n] << " " << nodes[3*n+1] << " " << nodes[3*n+2] << "\n";
    }
    fp << "  </DataArray>\n";
    fp << "</Points>\n";
    // ----------------------------------------------------------------

    // CELLS ----------------------------------------------------------
    fp << "<Cells>\n";
    // CONNECTIVITY
    fp << "  <DataArray type=\"" << Cell_conn_t_description << "\" Name=\"connectivity\" format=\"ascii\">\n";
    
    for (Cell_offs_t c = 0; c < n_cells; ++c)
    {
        const Cell_offs_t n_comp = cell_offset[c+1]-cell_offset[c];
        fp << "  ";
        for (Cell_offs_t k = 0; k < (n_comp-1); ++k)
        {
            fp << cell_conn[cell_offset[c]+k] << " ";
        }
        fp << cell_conn[cell_offset[c]+n_comp-1] << "\n";
    }
    fp << "  </DataArray>\n";

    // OFFSET
    fp << "  <DataArray type=\"" << Cell_offs_t_description << "\" Name=\"offsets\" format=\"ascii\">\n";
    fp << "  ";
    {
        Cell_offs_t c = 1;
        while (c < (n_cells+1))
        {
            fp << cell_offset[c];
            c += 1;
            fp << ((c%cell_offset_stride == 0)? "\n  " : " ");
        }
    }
    fp << "\n  </DataArray>\n";

    // TYPE
    fp << "  <DataArray type=\"" << Cell_type_t_description << "\" Name=\"types\" format=\"ascii\">\n";
    fp << "  ";
    {
        Cell_offs_t c = 0;
        while (c < n_cells)
        {
            fp << (int) cell_type[c];
            c += 1;
            fp << ((c%cell_type_stride == 0)? "\n  " : " ");
        }
    }
    fp << "\n  </DataArray>\n";

    fp << "</Cells>\n";
    // ----------------------------------------------------------------

    // CELLS DATA -----------------------------------------------------
    fp << "<CellData Scalars=\"Scalars\">\n";
    for (int f = 0; f < n_cell_fields; ++f)
    {
        fp << "  <DataArray type=\"" << Int_t_description << "\" Name=\""+cell_field_name[f]+"\" format=\"ascii\">\n";
        fp << "  ";
        {
            Cell_offs_t c = 0;
            while (c < n_cells)
            {
                fp << cell_field[f][c];
                c += 1;
                fp << ((c%cell_field_stride == 0)? "\n  " : " ");
            }
        }
        fp << "\n  </DataArray>\n";
    }
    fp << "</CellData>" <<"\n";
    // ----------------------------------------------------------------

    // POINTS DATA ----------------------------------------------------
    fp << "<PointData Scalars=\"Scalars\">\n";
    for (int f = 0; f < n_nodal_fields; ++f)
    {
        fp << "  <DataArray type=\"" << Float_t_description << "\" Name=\""+nodal_field_name[f]+"\" format=\"ascii\">\n";
        fp << "  ";
        {
            Cell_conn_t n = 0;
            while (n < n_nodes)
            {
                fp << nodal_field[f][n];
                n += 1;
                fp << ((n%nodal_field_stride == 0)? "\n  " : " ");
            }
        }
        fp << "\n  </DataArray>\n";
    }
    fp << "</PointData>" <<"\n";
    // ----------------------------------------------------------------

    // CLOSING -------------------
    fp << "</Piece>\n";
    fp << "</UnstructuredGrid>\n";
    fp << "</VTKFile>\n";
    // ---------------------------
}

void print_unstructured_data_binary(std::ofstream & fp,
                                    const Cell_conn_t & n_nodes,
                                    const Cell_offs_t & n_cells,
                                    const std::vector<Float_t> & nodes,
                                    const std::vector<Cell_conn_t> & cell_conn,
                                    const std::vector<Cell_offs_t> & cell_offset,
                                    const std::vector<Cell_type_t> & cell_type,
                                    const std::vector<std::vector<Int_t>> & cell_field,
                                    const std::vector<std::string> & cell_field_name,
                                    const std::vector<std::vector<Float_t>> & nodal_field,
                                    const std::vector<std::string> & nodal_field_name)
{

}

}
}