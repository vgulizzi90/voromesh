// utils_voropp.cpp

#include "utils_voropp.hpp"

namespace voromesh
{
namespace voro_utils
{

inline int cycle_up(int a,int nup) {return a==nup-1?0:a+1;}
inline int cycle_down(int a,int nup) {return a==0?nup-1:a-1;}

bool find_face(voro::voronoicell_neighbor &c,int &i,int &j,int q) {
    for(i=0;i<c.p;i++) for(j=0;j<c.nu[i];j++) if(c.ne[i][j]==q) return false;
    return true;
}

void face_neighbors(voro::voronoicell_neighbor &c,int &i,int &j,int q,std::vector<int> &vi) {
    int *nu=c.nu,**ed=c.ed,**ne=c.ne,k,l,m;
    vi.clear();
    vi.push_back(ne[i][cycle_down(j,nu[i])]);
    k=ed[i][j];
    l=cycle_up(ed[i][nu[i]+j],nu[k]);
    do {
        vi.push_back(ne[k][cycle_down(l,nu[k])]);
        m=ed[k][l];
        l=cycle_up(ed[k][nu[k]+l],nu[m]);
        k=m;
    } while(k!=i);
}

double sc_prod(voro::voronoicell_neighbor &c,double x,double y,double z,double rsq,int n,double &ans) {
    double *pp=c.pts+n+(n<<1);
    ans=*(pp++)*x;
    ans+=*(pp++)*y;
    ans+=*(pp++)*z-rsq;
    return ans;
}

/** Splits a face by adding in an additional edge across it.
 * \param[in] c the Voronoi cell to split.
 * \param[in] q the neighbor ID of the face to split.
 * \param[in] q2 the new neighbor ID of the additional face.
 * \param[in] (x,y,z) the normal vector of the cutting plane to split by.
 * \param[in] rsq the distance along the normal vector of the cutting plane. */
void split_face(voro::voronoicell_neighbor &c,int q,int q2,double x,double y,double z,double rsq) {
    int i,j=0,jj,k,**ed=c.ed,*nu=c.nu,**ne=c.ne,oi,oj,&p=c.p;
    double *pts=c.pts,ans,oans,r,ro;

    // Find an edge that goes along the side of the face to consider
    find_face(c,i,j,q);
    int i_=i,j_=j;

    // Cycle around the edges of the face until finding a vertex that is
    // inside the half space defined by the plane.
    while(sc_prod(c,x,y,z,rsq,i,ans)>0) {
        k=ed[i][j];
        j=cycle_up(ed[i][nu[i]+j],nu[k]);
        if(k==i_&&j==j_) {
            fputs("Error 1\n",stderr);
            exit(1);
        }
        i=k;
    }

    // Cycle around the edges of the face until finding a vertex that is
    // outside the half space defined by the plane.
    i_=i;j_=j;
    do {
        oi=i;oj=j;oans=ans;
        i=ed[oi][oj];
        j=cycle_up(ed[oi][nu[oi]+oj],nu[i]);
        if(i==i_&&j==j_) {
            fputs("Error 2\n",stderr);
            exit(1);
        }
    } while(sc_prod(c,x,y,z,rsq,i,ans)<=0);

    // Allocate memory for two new vertices
    if(p==c.current_vertices) c.add_memory_vertices(c);
    if(c.mec[3]==c.mem[3]) c.add_memory(c,3,NULL);
    ed[p]=c.mep[3]+7*c.mec[3];
    ne[p]=c.mne[3]+3*c.mec[3]++;
    ed[p+1]=c.mep[3]+7*c.mec[3];
    ne[p+1]=c.mne[3]+3*c.mec[3]++;

    // We now have an intersected edge. Create a new vertex along the edge.
    r=oans/(oans-ans);ro=1-r;
    pts[3*p]=pts[3*i]*r+pts[3*oi]*ro;
    pts[3*p+1]=pts[3*i+1]*r+pts[3*oi+1]*ro;
    pts[3*p+2]=pts[3*i+2]*r+pts[3*oi+2]*ro;

    // Set the edge and neighbor information for the new vertex, and update
    // the relations with neighboring vertices
    jj=ed[oi][oj+nu[oi]];
    nu[p]=3;
    ed[p][0]=p+1;
    ed[p][1]=i;
    ed[p][2]=oi;
    ed[p][3]=0;
    ed[p][4]=jj;
    ed[p][5]=oj;
    ed[p][6]=p;
    ed[oi][oj]=p;
    ed[oi][oj+nu[oi]]=2;
    ed[i][jj]=p;
    ed[i][jj+nu[i]]=1;
    ne[p][0]=q;
    ne[p][1]=q2;
    ne[p++][2]=ne[oi][cycle_up(oj,nu[oi])];

    // Cycle around all of the vertices that lie outside of the half space.
    // In addition, switch the neighbor information to the new neighbor ID.
    i_=i;j_=j;
    do {
        ne[i][j]=q2;
        oi=i;oj=j;oans=ans;
        i=ed[oi][oj];
        j=cycle_up(ed[oi][nu[oi]+oj],nu[i]);
        if(i==i_&&j==j_) {
            fputs("Error 3\n",stderr);
            exit(1);
        }
    } while(sc_prod(c,x,y,z,rsq,i,ans)>0);

    // We now have another intersected edge. Create a new vertex.
    r=oans/(oans-ans);ro=1-r;
    pts[3*p]=pts[3*i]*r+pts[3*oi]*ro;
    pts[3*p+1]=pts[3*i+1]*r+pts[3*oi+1]*ro;
    pts[3*p+2]=pts[3*i+2]*r+pts[3*oi+2]*ro;

    // Set the edge and neighbor information for the new vertex, and update
    // the relations with neighboring vertices
    nu[p]=3;
    jj=ed[oi][oj+nu[oi]];
    ed[p][0]=p-1;
    ed[p][1]=i;
    ed[p][2]=oi;
    ed[p][3]=0;
    ed[p][4]=jj;
    ed[p][5]=oj;
    ed[p][6]=p;
    ed[oi][oj]=p;
    ed[oi][oj+nu[oi]]=2;
    ed[i][jj]=p;
    ed[i][jj+nu[i]]=1;
    ne[p][0]=q2;
    ne[p][1]=q;
    ne[p++][2]=ne[oi][cycle_up(oj,nu[oi])];
}

}
}