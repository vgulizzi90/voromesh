// utils_math.hpp

#ifndef VOROMESH_MATH_HPP_
#define VOROMESH_MATH_HPP_

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

namespace voromesh
{
namespace math
{

#define THIRD (1.0/3.0)
#define SIXTH (1.0/6.0)

inline
double L2_norm(const double * v, const size_t n)
{
    double tmp;
    tmp = 0.0;
    for (size_t i = 0; i < n; ++i) tmp += v[i]*v[i];

    return std::sqrt(tmp);
}

inline
double L2_distance(const double * v1, const double * v2, const size_t n)
{
    double tmp;
    tmp = 0.0;
    for (size_t i = 0; i < n; ++i) tmp += (v1[i]-v2[i])*(v1[i]-v2[i]);

    return std::sqrt(tmp);
}

// res = v1 x v2
inline
void cross(const double * v1, const double * v2, double * res)
{
    res[0] = v1[1]*v2[2]-v1[2]*v2[1];
    res[1] = v1[2]*v2[0]-v1[0]*v2[2];
    res[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

// res = (A-O) x (B-O)
inline
void cross(const double * O, const double * A, const double * B, double * res)
{
    res[0] = -O[2]*A[1]+O[1]*A[2]+O[2]*B[1]-A[2]*B[1]-O[1]*B[2]+A[1]*B[2];
    res[1] = +O[2]*A[0]-O[0]*A[2]-O[2]*B[0]+A[2]*B[0]+O[0]*B[2]-A[0]*B[2];
    res[2] = -O[1]*A[0]+O[0]*A[1]+O[1]*B[0]-A[1]*B[0]-O[0]*B[1]+A[0]*B[1];
}

inline
void tri_centroid_2d(const double * v1, const double * v2, const double * v3, double * G)
{
    G[0] = THIRD*(v1[0]+v2[0]+v3[0]);
    G[1] = THIRD*(v1[1]+v2[1]+v3[1]);
}

inline
void tri_centroid(const double * v1, const double * v2, const double * v3, double * G)
{
    G[0] = THIRD*(v1[0]+v2[0]+v3[0]);
    G[1] = THIRD*(v1[1]+v2[1]+v3[1]);
    G[2] = THIRD*(v1[2]+v2[2]+v3[2]);
}

inline
double tri_area_3d(const double * v1, const double * v2, const double * v3)
{
    double tmp_v[3], tmp;

    cross(v1, v2, v3, tmp_v);
    tmp = std::sqrt(tmp_v[0]*tmp_v[0]+tmp_v[1]*tmp_v[1]+tmp_v[2]*tmp_v[2]);
    tmp = 0.5*tmp;

    return tmp;
}

inline
void tet_centroid(const double * v1, const double * v2, const double * v3, const double * v4, double * G)
{
    G[0] = 0.25*(v1[0]+v2[0]+v3[0]+v4[0]);
    G[1] = 0.25*(v1[1]+v2[1]+v3[1]+v4[1]);
    G[2] = 0.25*(v1[2]+v2[2]+v3[2]+v4[2]);
}

inline
double tet_volume(const double * v1, const double * v2, const double * v3, const double * v4)
{
    double tmp_v[3], tmp;

    cross(v1, v2, v3, tmp_v);
    tmp = (v4[0]-v1[0])*tmp_v[0]+(v4[1]-v1[1])*tmp_v[1]+(v4[2]-v1[2])*tmp_v[2];
    tmp = SIXTH*std::abs(tmp);

    return tmp;
}

inline
void linspace(const double xlo, const double xhi, const int N, double * S)
{
    const double dx = (xhi-xlo)/(N-1);
    S[0] = xlo;
    S[N-1] = xhi;
    for (int n = 1; n < (N-1); ++n)
    {
        S[n] = xlo+n*dx;
    }
}

// DELAUNAY TRIANGULATION #############################################
void delaunay_2d(std::vector<double> & points, std::vector<int> & tri, std::vector<int> & tri_edges,
                 const double min_quality = 1.0e-5);
void delaunay_3d(std::vector<double> & points, std::vector<int> & tet, std::vector<int> & tet_edges,
                 const double min_quality = 1.0e-5);
// ####################################################################

// MATRIX OPERATION ###################################################
namespace matrix
{

void print_inline(const int N, const double * S);
void print_inline_int(const int N, const int * S);

inline
void transpose(const int Nr, const int Nc, const double * S, double * D)
{
    for (int c = 0; c < Nc; ++c)
    for (int r = 0; r < Nr; ++r)
        D[c+r*Nc] = S[r+c*Nr];
}

// D <- S1*S2
inline
void multiply(const int Nr, const int Nc, const int Nrhs, const double * S1, const double * S2, double * D)
{
    int N = Nr*Nrhs;
    for (int n = 0; n < N; ++n) D[n] = 0.0;

    for (int k = 0; k < Nrhs; ++k)
    for (int r = 0; r < Nr; ++r)
    for (int c = 0; c < Nc; ++c)
        D[r+k*Nr] += S1[r+c*Nr]*S2[c+k*Nc];
}

}
// ####################################################################


}
}

#endif