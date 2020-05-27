// utils_voropp.hpp

#ifndef VOROMESH_UTILS_VOROPP_HPP_
#define VOROMESH_UTILS_VOROPP_HPP_

#include <voro++.hh>

namespace voromesh
{
namespace voro_utils
{

bool find_face(voro::voronoicell_neighbor &c,int &i,int &j,int k);

void face_neighbors(voro::voronoicell_neighbor &c,int &i,int &j,int k,std::vector<int> &vi);

double sc_prod(voro::voronoicell_neighbor &c,double x,double y,double z,double rsq,int n,double &ans);

void split_face(voro::voronoicell_neighbor &c,int q,int q2,double x,double y,double z,double rsq);

}
}

#endif