// utils_vtk.hpp

#ifndef VOROMESH_VTK_HPP_
#define VOROMESH_VTK_HPP_

#include <vector>
#include <string>
#include <fstream>

#include "utils_io.hpp"

#define VTK_VERTEX 1
#define VTK_LINE 3
#define VTK_TRIANGLE 5
#define VTK_QUAD 9
#define VTK_TETRA 10

namespace voromesh
{
namespace vtk
{

typedef uint32_t Header_t;
typedef double Float_t;
typedef int Int_t;
typedef uint32_t Cell_conn_t;
typedef uint32_t Cell_offs_t;
typedef uint8_t Cell_type_t;

const std::string Header_t_description = "UInt32";
const std::string Float_t_description = "Float"+std::to_string(sizeof(Float_t)*8);
const std::string Int_t_description = "Int"+std::to_string(sizeof(Int_t)*8);
const std::string Cell_conn_t_description = "Int"+std::to_string(sizeof(Cell_conn_t)*8);
const std::string Cell_offs_t_description = "Int"+std::to_string(sizeof(Cell_offs_t)*8);
const std::string Cell_type_t_description = "UInt"+std::to_string(sizeof(Cell_type_t)*8);

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
                             const std::string & fmt = "ascii");
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
                                   const std::vector<std::string> & nodal_field_name);
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
                                    const std::vector<std::string> & nodal_field_name);

}
}

#endif