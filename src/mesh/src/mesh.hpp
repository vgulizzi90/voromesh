// mesh.hpp

#ifndef VOROMESH_MESH_HPP_
#define VOROMESH_MESH_HPP_

#include <vector>
#include <voro++.hh>

namespace voromesh
{

// MESH ###############################################################
struct Mesh
{
    // DATA MEMBERS ===================================================
    std::vector<double> nodes;
    std::vector<int> etype;
    std::vector<int> conn;
    // ================================================================
};
// ####################################################################

}

#endif