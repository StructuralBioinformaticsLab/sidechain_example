#ifndef _MWIS_H
#define _MWIS_H

#include "graph.h"

#ifdef __cplusplus
extern "C" {
#endif

void mwis(struct graph_t *graph , double *glpk_x, unsigned char *x);
void gradproj_CC (struct graph_t *graph , double * x);
void myestimation(struct graph_t *graph, double *lp_x, unsigned char *x);

void gradproj_singh (struct graph_t *graph , double * x);



#ifdef __cplusplus
}
#endif

#endif /* _MWIS_H */

