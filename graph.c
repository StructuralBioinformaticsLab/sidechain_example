/*
Copyright (c) 2014, Mohammad Moghadasi
Division of System Engineering, Boston University
Structural Bioinformatics Laboratory, Boston University
All rights reserved.

E-mail: <mohamad@bu.edu>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
- Neither the name of the author nor the names of its contributors may be used
  to endorse or promote products derived from this software without specific
  prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <libgen.h>
#include <stdarg.h>
#include <stdint.h>
#include _MOL_INCLUDE_
#include "main.h"
#include "parameters.h"
#include "rotamer_placement.h"
#include "pack.h"
#include "mwis.h"
#include "graph.h"

//________________________________________________________________________

// graph.c

//________________________________________________________________________

static void hbondeng_list( struct atomgrp *ag, double *energy, int *list1, int numOfAtomsList1, int *list2, int numOfAtomsList2 , double rc)
{

    double rc2 = rc * rc;

    ( *energy ) = 0;

    for(int i=0;i<numOfAtomsList1;i++) {
	int ai=list1[i];
        mol_atom *atom_i = &( ag->atoms[ ai ] );
        if ( !( atom_i->hprop & DONATABLE_HYDROGEN ) && !( atom_i->hprop & HBOND_ACCEPTOR ) ) continue;

        for( int j = 0; j < numOfAtomsList2; j++ ) {
            int aj = list2[j];
            mol_atom *atom_j = &( ag->atoms[ aj ] );

            if ( atom_i->hprop & DONATABLE_HYDROGEN )
              {
                if ( !( atom_j->hprop & HBOND_ACCEPTOR ) ) continue;

                ( *energy ) += get_pairwise_hbondeng_nblist( ag->atoms, ai, ag->atoms, aj, NULL, rc2, 1, 1.0 );
              }
            else
              {
                if ( !( atom_j->hprop & DONATABLE_HYDROGEN ) ) continue;

                ( *energy ) += get_pairwise_hbondeng_nblist( ag->atoms, aj, ag->atoms, ai, NULL, rc2, 1, 1.0 );
              }
          }
      }
}

static void get_list_noRes_atoms( struct atomgrp *ag , int resi_num , int* noRes_list , int* natoms_noRes)
{
        int nat;
        int first_at;
        int ind = 0 ;

        for( int res = 0 ; res < ag->nres ; res++){
		if( res != resi_num ){
			nat = ag->iares[res+1] - ag->iares[res];
			first_at = ag->iares[res];
			for( int k = 0  ; k < nat ; k++ )
                                noRes_list[ ind + k ] = first_at + k ;
			ind += nat ;
		}
                
        }

        natoms_noRes[0] = ind ;
}

static void get_list_noIF_atoms( struct atomgrp *ag , struct ifres_list* res_list , int* noIF_list , int* natoms_noIF)
{
	int nresi = ag->nres; 
	int nIFres = res_list->num_of_ifres;
	int nat;
	int indicator;
	int first_at;
	// int last_at;
	int ind = 0 ;

	for( int res = 0 ; res < nresi ; res++){

	        indicator = 1 ; // noIF res
		for(int j = 0 ; j < nIFres ; j++){
			if( res == res_list->ifres_num[j]){ 
				indicator = 0;
				j = nIFres ;
			}
		}
		if( indicator == 1){
			nat = ag->iares[res+1] - ag->iares[res];
			first_at = ag->iares[res]  ;	
			for( int k = 0  ; k < nat ; k++ )
				noIF_list[ ind + k ] = first_at + k ;
			ind += nat ;
		}
	}
	natoms_noIF[0] = ind ;

}

void graph_t_malloc( struct graph_t** graph , int MAX_NODE , int N_CLIQUES)
{
        *graph = (struct graph_t*) _mol_malloc (sizeof (struct graph_t));
        (*graph)->weight = (double*) _mol_malloc(MAX_NODE * sizeof(double));    // NEW
        (*graph)->adjacency_matrix = (unsigned char*) _mol_malloc(MAX_NODE * MAX_NODE * sizeof(unsigned char)); // New

	(*graph)->clique_list = (int**) _mol_malloc(N_CLIQUES * sizeof(int*));
	for(int i = 0 ; i < N_CLIQUES ; i++)
		(*graph)->clique_list[i] = (int*) _mol_malloc(MAX_NODE * sizeof(int));

	(*graph)->clique_size = (int*) _mol_malloc(N_CLIQUES * sizeof(int));

	(*graph)->cliquesof_i = (int**) _mol_malloc(MAX_NODE * sizeof(int*));
	for(int i = 0 ; i < MAX_NODE ; i++)
		(*graph)->cliquesof_i[i] = (int*) _mol_malloc(N_CLIQUES * sizeof(int));

	(*graph)->node_nClq = (int*) _mol_malloc(MAX_NODE * sizeof(int));
}

void weight_edit(struct graph_t* graph)
{
        // Re-ordering the energies

        LARGE = graph->weight[0];
        MINIM = graph->weight[0];

        for ( int i = 0 ; i < graph->vertex_count ; i++){
                LARGE = (graph->weight[i]>LARGE) ? graph->weight[i] : LARGE;
                MINIM = (graph->weight[i]<MINIM) ? graph->weight[i] : MINIM;
        }

        LARGE = LARGE + 5;

//	printf("Large = %f Min = %f\n", LARGE, MINIM);

        for ( int i = 0 ; i < graph->vertex_count ; i++){
                graph->weight[i] = LARGE - graph->weight[i] ;

        }
}


void Free_graph_t( struct graph_t** graph , int MAX_NODE, int N_CLIQUES )
{
	free((*graph)->adjacency_matrix);
	free((*graph)->weight);
	free(*graph);

	for(int i = 0 ; i < N_CLIQUES ; i++)
		free((*graph)->clique_list[i]);
	free((*graph)->clique_list);
	free((*graph)->clique_size);

	for(int i = 0 ; i < MAX_NODE ; i++)
		free((*graph)->cliquesof_i[i]);
	free((*graph)->cliquesof_i);
	free((*graph)->node_nClq);
}

//


void eng_linklist_gen( struct atomgrp *RecLigAg , struct ifres_list* res_list , struct node_energy* energy , struct link_list* ln_list , 			  struct rot_info* rotinf[])
{
        int resi; // residue id for i
        int resj; // residue id for j
        int list1[20];
        int list2[20];
        int list1_size;
        int list2_size;

        struct agsetup* ags;
        ags=malloc(sizeof(struct agsetup));

//      init_nblst(RecLigAg , ags);		
//      update_nblst(RecLigAg , ags);		

        // List of non-Interface atoms in the backbone for self energy purposes
        int noIF_list[RecLigAg->natoms];
        int natoms_noIF;
        get_list_noIF_atoms ( RecLigAg , res_list , noIF_list , &natoms_noIF);
        int noRes_list[RecLigAg->natoms];
        int natoms_noRes;

        // Initializaing the Graph Structs
//	#pragma omp parallel for
        for( int i = 0 ; i < res_list->num_of_ifres ; i++ ){
               for(int j = 0 ; j < res_list->num_of_ifres ; j++){

                        for( int ri = 0 ; ri < rotinf[i]->nrotamers + 1 ; ri++){
                                for( int rj = 0 ; rj < rotinf[j]->nrotamers + 1 ; rj++){

                                        energy->residues[i][j].rotamers[ri][rj]    = 0 ;
                                        ln_list->residues[i][j].rotamers[ri][rj]   = 0 ;
                                }
                        }
                }
        }

        energy->nresidue  = res_list->num_of_ifres ;
        ln_list->nresidue = res_list->num_of_ifres ;

        // Filling the Graph Structs
        double noIF_eng;
        double inter_eng;
        double hbond_eng;
        double rc = 5;
        double prob_eng;
        double prob1 ;

        for( int i = 0 ; i < res_list->num_of_ifres ; i++ ){
                energy->nrotamers[i]  = rotinf[i]->nrotamers ;
                ln_list->nrotamers[i] = rotinf[i]->nrotamers ;
                resi = res_list->ifres_num[i];
                get_list( RecLigAg, resi , list1 , &list1_size);
                get_list_noRes_atoms ( RecLigAg , resi , noRes_list , &natoms_noRes);

                for( int ri = 0 ; ri < rotinf[i]->nrotamers + 1 ; ri++){
                        apply_rot_to_atmgrp( RecLigAg, rotinf[i], resi , ri );  // +1 : since rotamer indices starts at 1 IMP ***

                        vdweng_inter( RecLigAg, rc, list1, list1_size, noIF_list, natoms_noIF, &noIF_eng);
                        prob1 = rotinf[i]->rot[ri].probability;
                        prob_eng = -1 * KT * log(prob1/P0);

                        energy->residues[i][i].rotamers[ri][ri]= vdw_noIF_coef*noIF_eng + prob1_coef*prob_eng; //+vdw_self_coef*self_eng;

                        if( (energy->residues[i][i].rotamers[ri][ri] < big_eng) || (ri == 0) ){
                           ln_list->residues[i][i].rotamers[ri][ri] = 1 ;

                           for(int j = i+1 ; j < res_list->num_of_ifres ; j++){
                                resj = res_list->ifres_num[j];
                                get_list( RecLigAg, resj , list2 , &list2_size);

                                for( int rj = 0 ; rj < rotinf[j]->nrotamers + 1 ; rj++){
                                        apply_rot_to_atmgrp( RecLigAg, rotinf[j], resj , rj );// +1 since rot indices starts at 1

                                        vdweng_inter( RecLigAg, rc, list1, list1_size, list2, list2_size, &inter_eng);
                                        hbond_eng = 0 ;
                                        hbondeng_list( RecLigAg , &hbond_eng, list1, list1_size, list2, list2_size , rc);
                                        energy->residues[i][j].rotamers[ri][rj]= vdw_inter_coef * inter_eng + hbond_coef * hbond_eng;
                                        energy->residues[j][i].rotamers[rj][ri]= energy->residues[i][j].rotamers[ri][rj];

                                        if(  (energy->residues[i][j].rotamers[ri][rj] < big_eng) || (ri == 0) || (rj == 0) ) {
                                                ln_list->residues[i][j].rotamers[ri][rj] = 1 ;
                                                ln_list->residues[j][i].rotamers[rj][ri] = 1 ;
                                        }
                              }
                                apply_rot_to_atmgrp( RecLigAg, rotinf[j], resj , 0 );
                          }
                      }
                 }
                apply_rot_to_atmgrp( RecLigAg, rotinf[i], resi , 0 );
        }
	free(ags);
}

//================================

void node_counter( struct ifres_list* res_list , struct link_list* ln_list , struct rot_info* rotinf[], int* nnodes)
{
   int node_cnt = 0 ;

   if( res_list->num_of_ifres == 2 ){        // Two IF residues : New Version of Graph :: Triple to One

        for( int i = 0 ; i < res_list->num_of_ifres ; i++ ){
                 for(int j = i+1 ; j < res_list->num_of_ifres ; j++){   // New Graph
                        for( int ri = 0 ; ri < rotinf[i]->nrotamers+1 ; ri++){
                                for( int rj = 0 ; rj < rotinf[j]->nrotamers+1 ; rj++){
                                        if( ln_list->residues[i][j].rotamers[ri][rj] == 1  && ln_list->residues[i][i].rotamers[ri][ri] == 1                                            && ln_list->residues[j][j].rotamers[rj][rj] == 1){
                                                node_cnt++;
                                        }
                                }
                        }
                }
        }

    }

    else {                                    // More than Two IF residues : Old Version of Graph
    	for( int i = 0 ; i < res_list->num_of_ifres ; i++ ){
                 for(int j = i ; j < res_list->num_of_ifres ; j++){     // Old Graph
                        for( int ri = 0 ; ri < rotinf[i]->nrotamers+1 ; ri++){
                                for( int rj = 0 ; rj < rotinf[j]->nrotamers+1 ; rj++){
                                        if( ln_list->residues[i][j].rotamers[ri][rj] == 1  && ln_list->residues[i][i].rotamers[ri][ri] == 1                                            && ln_list->residues[j][j].rotamers[rj][rj] == 1){
                                                node_cnt++;
                                        }
                                }
                        }
                }
        }
    }
    *nnodes = node_cnt;
}

//================================

void graph_construction( struct graph_t* graph , struct ifres_list* res_list , struct link_list* ln_list ,  struct node_energy* energy , 			    struct rot_info* rotinf[] ,	int node_info[][4] , int nnodes)
{
     graph->vertex_count = nnodes ;
     graph->ncliques = res_list->num_of_ifres + res_list->num_of_ifres*(res_list->num_of_ifres - 1)/2;

     int node_ind = 0 ;
     int clique_ind = -1 ;
     int clqnod_ind[graph->ncliques];
     int nodescliq_ind[nnodes];		// nodescliq_ind[i] : index for the cliques of node i 


     if( res_list->num_of_ifres == 2 ){        // Two IF residues : New Version of Graph :: Triple to One

//	#pragma omp parallel for     

        for( int i = 0 ; i < res_list->num_of_ifres ; i++ ){
               for(int j = i+1 ; j < res_list->num_of_ifres ; j++){     // New Graph
                        for( int ri = 0 ; ri < rotinf[i]->nrotamers+1 ; ri++){
                                for( int rj = 0 ; rj < rotinf[j]->nrotamers+1 ; rj++){
                                        if( ln_list->residues[i][j].rotamers[ri][rj] == 1 && ln_list->residues[i][i].rotamers[ri][ri] == 1
                                         && ln_list->residues[j][j].rotamers[rj][rj] == 1){

                                                graph->weight[node_ind] =       energy->residues[i][i].rotamers[ri][ri]
                                                                              + energy->residues[j][j].rotamers[rj][rj]
                                                                              + energy->residues[i][j].rotamers[ri][rj];

                                                node_info[node_ind][0] = i;
                                                node_info[node_ind][1] = j;
                                                node_info[node_ind][2] = ri;
                                                node_info[node_ind][3] = rj;
                                                node_ind++;
                                        }
                                }
                        }
                }
        }
     }


     else {                                    // More than Two IF residues : Old Version of Graph
//           #pragma omp parallel for
     	     for(int i = 0 ; i < nnodes ; i++)  nodescliq_ind[i] = 0;     

//	#pragma omp parallel for     
        for( int i = 0 ; i < res_list->num_of_ifres ; i++ ){
               for(int j = i ; j < res_list->num_of_ifres ; j++){       // Old Graph
		        clique_ind++;
			clqnod_ind[clique_ind] = 0;
                        for( int ri = 0 ; ri < rotinf[i]->nrotamers+1 ; ri++){
                                for( int rj = 0 ; rj < rotinf[j]->nrotamers+1 ; rj++){
                                        if( ln_list->residues[i][j].rotamers[ri][rj] == 1 && ln_list->residues[i][i].rotamers[ri][ri] == 1
                                         && ln_list->residues[j][j].rotamers[rj][rj] == 1){

                                                graph->weight[node_ind] =       energy->residues[i][j].rotamers[ri][rj] ;

                                                node_info[node_ind][0] = i;
                                                node_info[node_ind][1] = j;
                                                node_info[node_ind][2] = ri;
                                                node_info[node_ind][3] = rj;

						// Adding Clique Info. of K-cliques when K > 2
						graph->clique_list[clique_ind][clqnod_ind[clique_ind]] = node_ind;		
						graph->cliquesof_i[node_ind][nodescliq_ind[node_ind]] = clique_ind;

						nodescliq_ind[node_ind]++;
                                                node_ind++;
						clqnod_ind[clique_ind]++;


                                        }
                                }
                        }
                }
        }

     

	for(int i = 0 ; i < graph->ncliques ; i++)
		graph->clique_size[i] = clqnod_ind[i];

	for(int i = 0 ; i < nnodes ; i++)
		graph->node_nClq[i] = nodescliq_ind[i];
    }


        // adjacency_matrix generation

//      #pragma omp parallel for
        for(unsigned int i = 0 ; i < graph->vertex_count ; i++){
                graph_edge(graph, i, i) = 0;
                for(unsigned int j = i+1 ; j < graph->vertex_count ; j++){

                        if(  ( (node_info[i][0] == node_info[j][0]) & (node_info[i][2] != node_info[j][2]) ) |
                             ( (node_info[i][1] == node_info[j][1]) & (node_info[i][3] != node_info[j][3]) ) |
                             ( (node_info[i][0] == node_info[j][1]) & (node_info[i][2] != node_info[j][3]) ) |
                             ( (node_info[i][1] == node_info[j][0]) & (node_info[i][3] != node_info[j][2]))) {

                                graph_edge(graph, i, j) = 1;
                                graph_edge(graph, j, i) = 1;
                        }
                        else{
                                graph_edge(graph, i, j) = 0;
                                graph_edge(graph, j, i) = 0;
                        }
                }
        }

}

//=====================================

void add_2cliques_to_graph( struct graph_t* graph)
{

	int clique_ind = graph->ncliques;

	for(int i = 0 ; i < graph->vertex_count ; i++)
		for(int j = i+1 ; j < graph->vertex_count ; j++)
			if(graph_edge(graph, i, j)){
				graph->clique_list[clique_ind][0] = i;
				graph->clique_list[clique_ind][1] = j;			
				graph->cliquesof_i[i][graph->node_nClq[i]] = clique_ind;
				graph->cliquesof_i[j][graph->node_nClq[j]] = clique_ind;

				graph->clique_size[clique_ind] = 2;

				clique_ind++;
				graph->node_nClq[i]=graph->node_nClq[i]+1;
				graph->node_nClq[j]=graph->node_nClq[j]+1;
			}

        graph->ncliques = clique_ind;
}
//=====================================

//________________________________________________________________________

// dee.c

//________________________________________________________________________

#ifndef minimum
        #define minimum( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif


void dee( struct node_energy* energy , struct link_list* ln_list )
{


	// Evaluating all the Best & Worst scenarios for all rotamers of all residues

	float  best[energy->nresidue][MAX_ROT];
        float worst[energy->nresidue][MAX_ROT];

	float MAX;
	float MIN;


	for( int r = 0 ; r < energy->nresidue ; r++ ){
		for( int i = 0 ; i <= energy->nrotamers[r] ; i++){    // r_i : residue r * rotamer i

			best[r][i]  = 0;
			worst[r][i] = 0;
			
			best[r][i]  += energy->residues[r][r].rotamers[i][i];
			worst[r][i] += energy->residues[r][r].rotamers[i][i]; 

			for ( int s = 0 ; s < energy->nresidue ; s++ ){
			   if( s != r){

				MAX = energy->residues[r][s].rotamers[i][0];
				MIN = energy->residues[r][s].rotamers[i][0];

				for( int j = 0 ; j <= energy->nrotamers[s] ; j++){    // s_j : residue s * rotamer j

					MAX = (energy->residues[r][s].rotamers[i][j]>MAX) ? energy->residues[r][s].rotamers[i][j]:MAX;
					MIN = (energy->residues[r][s].rotamers[i][j]<MIN) ? energy->residues[r][s].rotamers[i][j]:MIN;
				}
			   	
				best[r][i]  += MIN ;
				worst[r][i] += MAX ;
			   }
			}
		}
	}

	// Elimination Part

	for( int r = 0 ; r < energy->nresidue ; r++ ){
		for( int i1 = 0 ; i1 <= energy->nrotamers[r] ; i1++){
			for( int i2 = 0 ; i2 <= energy->nrotamers[r] ; i2++ ){

				if( (best[r][i1] > worst[r][i2]) && (i1 != 0) ) {  // rotamer i2 beats rotamer i1

						ln_list->residues[r][r].rotamers[i1][i1]=0;	
						break;			
				}
			}
		}
	}

	for( int r = 0 ; r < energy->nresidue ; r++ ){
                for( int i = 0 ; i <= energy->nrotamers[r] ; i++){
			if( ln_list->residues[r][r].rotamers[i][i] == 0){		

				for( int s = 0 ; s < energy->nresidue ; s++ ){
					for( int j = 0 ; j <= energy->nrotamers[s] ; j++){
		
						ln_list->residues[r][s].rotamers[i][j] = 0;
					}
				}
			}
		}
	}

}

// Memory Allocation node_enegy

void node_energy_malloc( struct node_energy** energy , int MAX_RES1 ){  // NEW

        *energy = (struct node_energy*) _mol_malloc (sizeof (struct node_energy ));
        (*energy)->residues = (struct eng_residue**) _mol_malloc ( MAX_RES1  * sizeof (struct eng_residue*));

        for(int i = 0 ; i < MAX_RES1 ; i++)
                (*energy)->residues[i] = (struct eng_residue*) _mol_malloc ( MAX_RES1  * sizeof (struct eng_residue));

        for(int i = 0 ; i < MAX_RES1 ; i++)
                for(int j = 0 ; j < MAX_RES1 ; j++)
                        (*energy)->residues[i][j].rotamers = (double**) _mol_malloc( MAX_ROT * sizeof(double*) );

        for(int i = 0 ; i < MAX_RES1 ; i++)
                for(int j = 0 ; j < MAX_RES1 ; j++)
                        for(int k = 0 ; k < MAX_ROT ; k++)
                                (*energy)->residues[i][j].rotamers[k] =  (double*) _mol_malloc( MAX_ROT * sizeof(double) );

        (*energy)->nrotamers = malloc ( MAX_RES1 * sizeof(int) );

}


void Free_node_energy( struct node_energy** energy , int MAX_RES1 )
{
        free( (*energy)->nrotamers );

        for(int i = 0 ; i < MAX_RES1 ; i++){
                for(int j = 0 ; j < MAX_RES1 ; j++){
                        for(int k = 0 ; k < MAX_ROT ; k++){
                                free( (*energy)->residues[i][j].rotamers[k] );
                        }
                }
        }

        for(int i = 0 ; i < MAX_RES1 ; i++){
                for(int j = 0 ; j < MAX_RES1 ; j++){
                        free( (*energy)->residues[i][j].rotamers);
                }
        }

        for(int i = 0 ; i < MAX_RES1 ; i++){
                free( (*energy)->residues[i] );
        }

        free((*energy)->residues);
        free(*energy );
}

// Memory Allocation link_list

void link_list_malloc( struct link_list** ln_list , int MAX_RES1 )      // NEW
{
        *ln_list = (struct link_list*) _mol_malloc (sizeof (struct link_list));

        (*ln_list)->residues = (struct rotamer_pair**) _mol_malloc ( MAX_RES1  * sizeof (struct rotamer_pair*));

        for(int i = 0 ; i < MAX_RES1 ; i++)
                (*ln_list)->residues[i] = (struct rotamer_pair*) _mol_malloc ( MAX_RES1  * sizeof (struct rotamer_pair));

        for(int i = 0 ; i < MAX_RES1 ; i++)
                for(int j = 0 ; j < MAX_RES1 ; j++)
                        (*ln_list)->residues[i][j].rotamers =  (int**) _mol_malloc (MAX_ROT * sizeof(int*) );

        for(int i = 0 ; i < MAX_RES1 ; i++)
                for(int j = 0 ; j < MAX_RES1 ; j++)
                        for(int k = 0 ; k < MAX_ROT ; k++)
                                (*ln_list)->residues[i][j].rotamers[k] = (int*) _mol_malloc (MAX_ROT * sizeof(int));

        (*ln_list)->nrotamers = malloc ( MAX_RES1 * sizeof(int) );

}

//

void Free_link_list ( struct link_list** ln_list , int MAX_RES1 )
{
        free((*ln_list)->nrotamers);

        for(int i = 0 ; i < MAX_RES1 ; i++){
                for(int j = 0 ; j < MAX_RES1 ; j++){
                        for(int k = 0 ; k < MAX_ROT ; k++){
                                free( (*ln_list)->residues[i][j].rotamers[k] );
                        }
                }
        }

        for(int i = 0 ; i < MAX_RES1 ; i++){
                for(int j = 0 ; j < MAX_RES1 ; j++){
                        free( (*ln_list)->residues[i][j].rotamers );
                }
        }

        for(int i = 0 ; i < MAX_RES1 ; i++){

                free((*ln_list)->residues[i]);
        }
        free((*ln_list)->residues);
        free(*ln_list);
}

// vdweng_list.c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define CUTOFF_LINE_SCALE          0.6
#define CUTOFF_REP_SCALE           0.89

void vdweng_inter(struct atomgrp *ag, double rc, int *list1, int numOfAtomsList1, int *list2, int numOfAtomsList2, double *eng){

        int i, i1, j, i2;
        double ei, ri, ej, rj, x1, y1, z1, dx, dy, dz, eij;
        double d2, Rd12, Rd6, Rr12, Rr6, dr6;
        double rij,rij2;
//      double dven, g;
        struct atom *a1, *a2;
        double rc2=rc*rc;
        double a,b,d;
        double cutoff_line_2, cutoff_rep_2;
        double cutoff_rep_scale_2;
        double cutoff_line_scale_2, cutoff_line_scale_6, cutoff_line_scale_min_6, cutoff_line_scale_min_12;
        double *eng_rep = malloc(sizeof(double));
	double *eng_atr = malloc(sizeof(double));
        *eng_rep = 0;
        *eng_atr = 0;
	(*eng) = 0 ;

        cutoff_rep_scale_2       = CUTOFF_REP_SCALE*CUTOFF_REP_SCALE;
        cutoff_line_scale_2      = CUTOFF_LINE_SCALE*CUTOFF_LINE_SCALE;
        cutoff_line_scale_6      = pow(CUTOFF_LINE_SCALE,6);
        cutoff_line_scale_min_6  = pow(CUTOFF_LINE_SCALE,-6);
        cutoff_line_scale_min_12 = pow(CUTOFF_LINE_SCALE,-12);

        for(i=0;i<numOfAtomsList1;i++){
           i1=list1[i];
           a1=&(ag->atoms[i1]);
           ei=a1->eps;
           ri=a1->rminh;
           x1=ag->atoms[i1].X;
           y1=ag->atoms[i1].Y;
           z1=ag->atoms[i1].Z;

           for(j=0;j<numOfAtomsList2;j++){
              i2=list2[j];
              a2=&(ag->atoms[i2]);
              rj=a2->rminh;
              ej=a2->eps;
              dx=x1-ag->atoms[i2].X;
              dy=y1-ag->atoms[i2].Y;
              dz=z1-ag->atoms[i2].Z;

              d2=dx*dx+dy*dy+dz*dz;

              if(d2<rc2){
                 eij  = ei*ej;
                 rij  = ri+rj;
                 rij2 = rij*rij;//artem modification for rij
                 cutoff_line_2 = cutoff_line_scale_2*rij2;
                 cutoff_rep_2  = cutoff_rep_scale_2*rij2;

                 if(d2<cutoff_line_2){

                    Rr6=rij2/rc2;
                    Rr6=Rr6*Rr6*Rr6;
                    Rr12=Rr6*Rr6;
                    d = sqrt(d2);

                    a = cutoff_line_scale_min_12 - 2*cutoff_line_scale_min_6 + Rr6*(4.0-2*cutoff_line_scale_6*Rr6) + Rr12*(2*cutoff_line_scale_6*Rr6-3.0);
                    b = 12.0*(-cutoff_line_scale_min_12 + cutoff_line_scale_min_6 + cutoff_line_scale_6*Rr6*(Rr12-Rr6))/(CUTOFF_LINE_SCALE*rij);

                    (*eng_rep) += eij*(a + (d - CUTOFF_LINE_SCALE*rij)*b);
                    //eng = eij*(a + (d - cutoff_line_scale*rij)*b);
                    //dven = -eij*b/d;
                    //printf("line eng at %f: %f, deriv: %f\n",d,*eng_rep,dven);

                 }else if(d2<cutoff_rep_2){

                    Rd6=rij2/d2;
                    Rd6=Rd6*Rd6*Rd6;
                    Rd12=Rd6*Rd6;

                    Rr6=rij2/rc2;
                    Rr6=Rr6*Rr6*Rr6;
                    Rr12=Rr6*Rr6;

                    dr6=d2/rc2;
                    dr6=dr6*dr6*dr6;

                    (*eng_rep) += eij*(Rd12 - 2*Rd6 + Rr6*(4.0-2*dr6) + Rr12*(2*dr6-3.0));
                    //eng = eij*(Rd12 - 2*Rd6 + Rr6*(4.0-2*dr6) + Rr12*(2*dr6-3.0));
                    //dven = -eij*12*(-Rd12+Rd6+dr6*(Rr12-Rr6))/d2;

                 }else{

                    Rd6=rij2/d2;
                    Rd6=Rd6*Rd6*Rd6;
                    Rd12=Rd6*Rd6;

                    Rr6=rij2/rc2;
                    Rr6=Rr6*Rr6*Rr6;
                    Rr12=Rr6*Rr6;

                    dr6=d2/rc2;
                    dr6=dr6*dr6*dr6;

                    (*eng_atr) += eij*(Rd12 - 2*Rd6 + Rr6*(4.0-2*dr6) + Rr12*(2*dr6-3.0));
                    //dven = -eij*12*(-Rd12+Rd6+dr6*(Rr12-Rr6))/d2;

                 }

                 /*
                 g=dven*dx;
                 (a1->GX)+=g;
                 (a2->GX)-=g;
                 g=dven*dy;
                 (a1->GY)+=g;
                 (a2->GY)-=g;
                 g=dven*dz;
                 (a1->GZ)+=g;
                 (a2->GZ)-=g;
                 */
              }
           }
        }

        *eng = *eng_rep + *eng_atr;
        free(eng_rep);
        free(eng_atr);

}



void get_list( struct atomgrp *ag , int resid , int* list, int* list_size)
{
        list_size[0] = ag->iares[resid+1] - ag->iares[resid];
        int first_at = ag->iares[resid]  ;
	// int last_at  = ag->iares[resid+1] - 1;

        for ( int i = 0 ; i < list_size[0] ; i++)
                list[i] = first_at + i ;

}




//





