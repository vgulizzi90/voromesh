/*****************************************************************************/
/*                                                                           */
/*  (tricall.c)                                                              */
/*                                                                           */
/*  Example program that demonstrates how to call Triangle.                  */
/*                                                                           */
/*  Accompanies Triangle Version 1.6                                         */
/*  July 19, 1996                                                            */
/*                                                                           */
/*  This file is placed in the public domain (but the file that it calls     */
/*  is still copyrighted!) by                                                */
/*  Jonathan Richard Shewchuk                                                */
/*  2360 Woolsey #H                                                          */
/*  Berkeley, California  94705-1927                                         */
/*  jrs@cs.berkeley.edu                                                      */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/*  report()   Print the input or output.                                    */
/*                                                                           */
/*****************************************************************************/

#include "tricall.h"

void report(io, markers, reporttriangles, reportneighbors, reportsegments,
            reportedges, reportnorms)
struct triangulateio *io;
int markers;
int reporttriangles;
int reportneighbors;
int reportsegments;
int reportedges;
int reportnorms;
{
    int i, j;

    for (i = 0; i < io->numberofpoints; i++) {
        printf("Point %4d:", i);
        for (j = 0; j < 2; j++) {
            printf("  %.6g", io->pointlist[i * 2 + j]);
        }
        if (io->numberofpointattributes > 0) {
            printf("   attributes");
        }
        for (j = 0; j < io->numberofpointattributes; j++) {
            printf("  %.6g",
                   io->pointattributelist[i * io->numberofpointattributes + j]);
        }
        if (markers) {
            printf("   marker %d\n", io->pointmarkerlist[i]);
        } else {
            printf("\n");
        }
    }
    printf("\n");

    if (reporttriangles || reportneighbors) {
        for (i = 0; i < io->numberoftriangles; i++) {
            if (reporttriangles) {
                printf("Triangle %4d points:", i);
                for (j = 0; j < io->numberofcorners; j++) {
                    printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);
                }
                if (io->numberoftriangleattributes > 0) {
                    printf("   attributes");
                }
                for (j = 0; j < io->numberoftriangleattributes; j++) {
                    printf("  %.6g", io->triangleattributelist[i *
                                                               io->numberoftriangleattributes + j]);
                }
                printf("\n");
            }
            if (reportneighbors) {
                printf("Triangle %4d neighbors:", i);
                for (j = 0; j < 3; j++) {
                    printf("  %4d", io->neighborlist[i * 3 + j]);
                }
                printf("\n");
            }
        }
        printf("\n");
    }

    if (reportsegments) {
        for (i = 0; i < io->numberofsegments; i++) {
            printf("Segment %4d points:", i);
            for (j = 0; j < 2; j++) {
                printf("  %4d", io->segmentlist[i * 2 + j]);
            }
            if (markers) {
                printf("   marker %d\n", io->segmentmarkerlist[i]);
            } else {
                printf("\n");
            }
        }
        printf("\n");
    }

    if (reportedges) {
        for (i = 0; i < io->numberofedges; i++) {
            printf("Edge %4d points:", i);
            for (j = 0; j < 2; j++) {
                printf("  %4d", io->edgelist[i * 2 + j]);
            }
            if (reportnorms && (io->edgelist[i * 2 + 1] == -1)) {
                for (j = 0; j < 2; j++) {
                    printf("  %.6g", io->normlist[i * 2 + j]);
                }
            }
            if (markers) {
                printf("   marker %d\n", io->edgemarkerlist[i]);
            } else {
                printf("\n");
            }
        }
        printf("\n");
    }
}

/*****************************************************************************/
/*                                                                           */
/*  main()   Create and refine a mesh.                                       */
/*                                                                           */
/*****************************************************************************/

int tricall(int nop,REAL **ipoints,int* onoed,int **oedges,int **oedgesmarker,int* onotri,int **otriangles)
{
    struct triangulateio in, out;
    unsigned int k;

    /* Define input points. */

    in.numberofpoints = nop;
    in.numberofpointattributes = 0;
    in.pointlist = (REAL *) malloc(in.numberofpoints*2*sizeof(REAL));
    for (k = 0; k < in.numberofpoints*2; ++k)
        in.pointlist[k] = (*ipoints)[k];

    in.pointattributelist = (REAL *) NULL;

    in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));
    for (k = 0; k < in.numberofpoints; ++k)
        in.pointmarkerlist[k] = 0;

    in.numberofsegments = 0;
    in.numberofholes = 0;
    in.numberofregions = 0;
    in.regionlist = (REAL *) NULL;

/*  printf("Input point set:\n\n");
    report(&in, 1, 0, 0, 0, 0, 0);*/

    /* Make necessary initializations so that Triangle can return a */
    /*   triangulation in `out'.  */

    out.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
    /* Not needed if -N switch used or number of point attributes is zero: */
    out.pointattributelist = (REAL *) NULL;
    out.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
    out.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
    /* Not needed if -E switch used or number of triangle attributes is zero: */
    out.triangleattributelist = (REAL *) NULL;
    out.neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
    /* Needed only if segments are output (-p or -c) and -P not used: */
    out.segmentlist = (int *) NULL;
    /* Needed only if segments are output (-p or -c) and -P and -B not used: */
    out.segmentmarkerlist = (int *) NULL;
    out.edgelist = (int *) NULL;             /* Needed only if -e switch used. */
    out.edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */

    /* Triangulate the points.  Switches are chosen to read and write a  */
    /*   PSLG (p), preserve the convex hull (c), number everything from  */
    /*   zero (z), assign a regional attribute to each element (A), and  */
    /*   produce an edge list (e), a Voronoi diagram (v), and a triangle */
    /*   neighbor list (n).                                              */

    triangulate("pczAenQ", &in, &out, (struct triangulateio *) NULL);

/*  printf("Initial triangulation:\n\n");
    report(&out, 1, 1, 1, 1, 1, 0);*/

    /* COPY RESULTS INTO oedges AND otriangles */
    (*onoed) = out.numberofedges;
    (*onotri) = out.numberoftriangles;
    (*oedges) = (int *) malloc(2*out.numberofedges*sizeof(int));
    (*oedgesmarker) = (int *) malloc(out.numberofedges*sizeof(int));
    (*otriangles) = (int *) malloc(3*out.numberoftriangles*sizeof(int));

    for (k = 0; k < out.numberofedges*2; ++k)
        (*oedges)[k] = out.edgelist[k];

    for (k = 0; k < out.numberofedges; ++k)
        (*oedgesmarker)[k] = out.edgemarkerlist[k];

    for (k = 0; k < out.numberoftriangles*3; ++k)
        (*otriangles)[k] = out.trianglelist[k];

    /* Free all allocated arrays, including those allocated by Triangle. */

    free(in.pointlist);
    free(in.pointattributelist);
    free(in.pointmarkerlist);
    free(in.regionlist);
    free(out.pointlist);
    free(out.pointattributelist);
    free(out.pointmarkerlist);
    free(out.trianglelist);
    free(out.triangleattributelist);
    free(out.neighborlist);
    free(out.segmentlist);
    free(out.segmentmarkerlist);
    free(out.edgelist);
    free(out.edgemarkerlist);

    return 0;
}
