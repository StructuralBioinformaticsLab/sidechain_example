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
#include <string.h>
#include <time.h>
#include <math.h>
#include <libgen.h>
#include <stdarg.h>
#include <unistd.h>
#include <stdint.h>
#include <limits.h>
#include _MOL_INCLUDE_
#include "main.h"
#include "parameters.h"
#include "rotamers.h"
#include "rotamer_placement.h"
#include "graph.h"
#include "mwis.h"
#include "omp.h"
#include "pack.h"
bool iprint = false;

static void reslist_clustering( struct atomgrp* RecLigAg , struct ifres_list* res_list , struct ifres_list* resList_C[] , int* nclusters)
{
        int i,j;
        int resi , resj;
        int* list1 = malloc( 20 * sizeof(int));
        int* list2 = malloc( 20 * sizeof(int));
        int list1_size;
        int list2_size;
        double rc = 5 ;
        double inter_eng;

	int matrix[res_list->num_of_ifres][res_list->num_of_ifres];

//	#pragma omp parallel for
	for(i = 0 ; i < res_list->num_of_ifres ; i++)
		for(j = 0 ; j < res_list->num_of_ifres ; j++)
			matrix[i][j]=0;

//	#pragma omp parallel for private (j, resi, resj)
	for(i = 0 ; i < res_list->num_of_ifres ; i++){
		
                resi = res_list->ifres_num[i];
                get_list( RecLigAg , resi , list1 , &list1_size);

		for(j = i+1 ; j < res_list->num_of_ifres ; j++){	       	                
			resj = res_list->ifres_num[j];
                        get_list( RecLigAg , resj , list2 , &list2_size);		

                        vdweng_inter( RecLigAg , rc, list1, list1_size, list2, list2_size, &inter_eng);

			if( inter_eng > 100){		// default 100
				matrix[i][j]=1;
				matrix[j][i]=1;
			}
		}
	}


        int check_list[res_list->num_of_ifres];                 // check_list[i] = 1 :: IFres[i] has been assigned to a cluster
	int cluster_atom[res_list->num_of_ifres];
	int assined_to_clus[res_list->num_of_ifres];
	int clus_size[res_list->num_of_ifres];
	int new_clus[res_list->num_of_ifres];

//	#pragma omp parallel for
	for(i = 0 ; i < res_list->num_of_ifres ; i++){
		check_list[i] = 0 ;
		cluster_atom[i]=0;
		assined_to_clus[i]=-1;
	        clus_size[i]=0;
		new_clus[i]=1;
	}

	int clus_ind = 0 ;
	int atom_ind = 0 ;

//	#pragma omp parallel for private(j)	// Segmentation Fault if Parallel

	for(i = 0 ; i < res_list->num_of_ifres ; i++){

		for(j = 0 ; j < res_list->num_of_ifres ; j++)
			if(matrix[i][j]==1)
				if(check_list[j]==1){	
					new_clus[i]=0;   // i should not start a cluster
                                        assined_to_clus[i] = assined_to_clus[j];
					break;
				}

		if((check_list[i] == 0) && (new_clus[i] == 1)){

			atom_ind = cluster_atom[clus_ind];
                        resList_C[clus_ind]->ifres_num[atom_ind] = res_list->ifres_num[i];
                        cluster_atom[clus_ind]++;
                        assined_to_clus[i]=clus_ind;
			clus_ind++;
                        check_list[i] = 1;	
			clus_size[assined_to_clus[i]]++;
			 
 			for(j = 0 ; j < res_list->num_of_ifres ; j++)
                        	if(matrix[i][j] == 1)
                                        if(check_list[j]==0){
	                                        resList_C[assined_to_clus[i]]->ifres_num[cluster_atom[assined_to_clus[i]]] = 																res_list->ifres_num[j];
                                                cluster_atom[assined_to_clus[i]]++;
                                                check_list[j] = 1;
						clus_size[assined_to_clus[i]]++;
						assined_to_clus[j] = assined_to_clus[i];
                                        }			
		}
	
		else if((check_list[i] == 0) && (new_clus[i] == 0)){

                        atom_ind = cluster_atom[assined_to_clus[i]];			
                        resList_C[assined_to_clus[i]]->ifres_num[atom_ind] = res_list->ifres_num[i];
                        cluster_atom[assined_to_clus[i]]++;
                        check_list[i] = 1;
			clus_size[assined_to_clus[i]]++;

                        for(j = 0 ; j < res_list->num_of_ifres ; j++)
                                if(matrix[i][j] == 1)
                                        if(check_list[j]==0){
                                                resList_C[assined_to_clus[i]]->ifres_num[cluster_atom[assined_to_clus[i]]] =                                                                                                                            res_list->ifres_num[j];
                                                cluster_atom[assined_to_clus[i]]++;
                                                check_list[j] = 1;
						clus_size[assined_to_clus[i]]++;
						assined_to_clus[j] = assined_to_clus[i];
                                        }		
		}				

		else if(check_list[i] == 1)	
			for( j = 0 ; j < res_list->num_of_ifres ; j++)
                                if(matrix[i][j] == 1)
                                        if(check_list[j]==0){
                                                resList_C[assined_to_clus[i]]->ifres_num[cluster_atom[assined_to_clus[i]]] =                                                                                                                            res_list->ifres_num[j];
                                                cluster_atom[assined_to_clus[i]]++;
                                                check_list[j] = 1;
						clus_size[assined_to_clus[i]]++;
						assined_to_clus[j] = assined_to_clus[i];
                                        }
	}


	nclusters[0] = clus_ind ;

//	#pragma omp parallel for		// + > -
	for(i = 0 ; i < clus_ind ; i++){
	        resList_C[i]->num_of_ifres = clus_size[i];
                array_sort( resList_C[i]->ifres_num , clus_size[i] );
	}

	free(list1);
	free(list2);	
}

static void reslist_nnodes( struct atomgrp* RecLigAg , struct ifres_list* res_list , int* nnodes, struct node_energy* energy                        , struct link_list* ln_list, struct rot_info* rotinf[] )
{
	// Generating Energy Graph & Link List

        eng_linklist_gen( RecLigAg , res_list , energy , ln_list , rotinf );    
        int nnodes0;
        node_counter( res_list , ln_list , rotinf, &nnodes0);   

	// Dead-End Elimination

        if( (DEE_i == 1) && (res_list->num_of_ifres > 2) )
                dee(energy , ln_list); 
	
        int nnodes1;
        node_counter( res_list , ln_list , rotinf, &nnodes1);
	*nnodes = nnodes1;
}

static void pack_reslist( struct atomgrp* RecLigAg , struct ifres_list* res_list , struct rot_info* rotinf_in, int* indexN)
{
	clock_t str, fin;

	str = clock();	////

        int i,j,k,l;

        int nnodes;
        int nresIF = res_list->num_of_ifres;

        struct node_energy* energy;
        struct link_list* ln_list;
        node_energy_malloc( &energy , nresIF );
        link_list_malloc( &ln_list, nresIF );

        struct rot_info* rotinf[nresIF];
        for(i = 0 ; i < nresIF ; i++)
                rot_info_malloc(&rotinf[i]);


        int i_total;
//	#pragma omp parallel for private(k,l)		// - > + 
        for(i = 0 ; i < nresIF ; i++){
                i_total = indexN[i];
                rotinf[i]->res_num = rotinf_in[i_total].res_num;
                rotinf[i]->natoms  = rotinf_in[i_total].natoms;
		rotinf[i]->nrotamers = rotinf_in[i_total].nrotamers;	
                for(k = 0 ; k <= rotinf[i]->nrotamers ; k++){
                        rotinf[i]->rot[k].probability = rotinf_in[i_total].rot[k].probability;
                        for(l = 0 ; l < rotinf[i]->natoms ; l++){
                                rotinf[i]->rot[k].X[l] = rotinf_in[i_total].rot[k].X[l];
                                rotinf[i]->rot[k].Y[l] = rotinf_in[i_total].rot[k].Y[l];
                                rotinf[i]->rot[k].Z[l] = rotinf_in[i_total].rot[k].Z[l];
                        }
                }
        }

	// Filter Out Low Probable Rotamers

	double sum_prob;

	#pragma omp parallel for private(k)
	for(i = 0 ; i < nresIF ; i++){
	   sum_prob = 0;
	   if(rotinf[i]->nrotamers > 1)
		for(k = 1 ; k <= rotinf[i]->nrotamers ; k++){
			sum_prob +=  rotinf[i]->rot[k].probability;
                        if(sum_prob >= Prob_Sum){

				rotinf[i]->nrotamers = k;
				break;
			}
		}
	}

	if(0)
           for(i = 0 ; i < nresIF ; i++)
               printf("nrotamers[%d]=%d\n",i,rotinf[i]->nrotamers);


	// Probability Correction

        double prob_sum;
        for(i = 0 ; i < res_list->num_of_ifres ; i++){

          if(rotinf[i]->nrotamers == 0)
		rotinf[i]->rot[0].probability = 1 ;

	  else if( UNB_i && rotinf[i]->nrotamers > 0 ){
               rotinf[i]->rot[0].probability = nb_rot_prob ;
               for(j = 1 ; j <= rotinf[i]->nrotamers ; j++)
                       prob_sum = prob_sum + rotinf[i]->rot[j].probability ;
               for(j = 1 ; j <= rotinf[i]->nrotamers ; j++)
                        rotinf[i]->rot[j].probability = rotinf[i]->rot[j].probability * (1 - nb_rot_prob ) / prob_sum ;
           }
        }

	reslist_nnodes( RecLigAg , res_list , &nnodes , energy , ln_list , rotinf );

	if(0)
           printf("Graph size after DEE = %d\n", nnodes );


        // Edit the Res_List for Extremely Large Clusters

	double rotamer_resuction_rate = .5;

        if( nnodes > Global_Max_Node ){
                while( nnodes > Global_Max_Node ){
			for(i = 0 ; i < res_list->num_of_ifres ; i++)
				rotinf[i]->nrotamers = floor( rotinf[i]->nrotamers * rotamer_resuction_rate) ;//+ 1;
			reslist_nnodes( RecLigAg , res_list, &nnodes , energy , ln_list , rotinf );
			if(0) printf("Graph size after Reduction = %d\n", nnodes );
		}													                
	}

	if(0) printf("Graph size = %d\n", nnodes );

	if(0)
	   for(i = 0 ; i < res_list->num_of_ifres ; i++)
	        printf("Final Number of Rotamers %s = %d\n", RecLigAg->idres[res_list->ifres_num[i]] ,rotinf[i]->nrotamers);	
	
	// To Generate the Cliques

	int ncliques_mal = nresIF + nresIF*(nresIF - 1)/2 + nnodes*(nnodes - 1)/2 ;




        fin = clock();
        if(iprint) printf("Part 1 : %f   sec\n", ((double)(fin-str)/CLOCKS_PER_SEC));



	str = clock();



	// Graph Generation Based on Res_List

        struct graph_t* graph ;
	        
	if(0) printf("nnodes = %d\n",nnodes);

        graph_t_malloc( &graph , nnodes, ncliques_mal );
        int node_info[nnodes][4];               // DICTIONARY : for node numbers conversion into rotamer data
        graph_construction( graph , res_list , ln_list , energy , rotinf , node_info , nnodes); // OK
        weight_edit(graph);             // Scale MAX - W(i) & Normalize W(i)'s to 1     // OK
        Free_link_list ( &ln_list , nresIF );
        Free_node_energy( &energy , nresIF );

	// Adding 2-cliques to the graph structure
	
	add_2cliques_to_graph(graph);


        fin = clock();
        if(iprint) printf("Part 2 : %f   sec\n", ((double)(fin-str)/CLOCKS_PER_SEC));



	str = clock();


        //  Solving the LP_Relaxation by Gradient Projection of Clique Constrained Problem

//      int MWIS_mode = 1;      // 0: my_MWIS   1: singh_MWIS

        double x_gp[graph->vertex_count];

        /*~~~~~ TWO MODES OF MWIS ALGORITHM ~~~~~*/

        clock_t start, finish;
        start = clock();

        if(!MWIS_mode)
           gradproj_CC (graph , x_gp);
        else
           gradproj_singh(graph, x_gp);

        finish = clock();
        if(iprint) printf("MWIS Time = %f\n",((double)(finish-start)/CLOCKS_PER_SEC));

        fin = clock();
        if(iprint) printf("Part 3.a : %f   sec\n", ((double)(fin-str)/CLOCKS_PER_SEC));

	str = clock();

        /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	                        


	// Estimation 

	unsigned char* x  = malloc( graph->vertex_count );	
	myestimation( graph, x_gp , x);

        // Finding the indices of the found nodes

        int node_list[nresIF];
        int ind;
	#pragma omp parallel for private(ind)
        for(i = 0 ; i < graph->vertex_count ; i++)
                if(x[i] == 1){
                      ind = i;
                      node_list[node_info[ind][0]] = node_info[ind][2]  ;
                      node_list[node_info[ind][1]] = node_info[ind][3]  ;
                }
	free(x);


        fin = clock();
        if(iprint) printf("Part 3.b : %f   sec\n", ((double)(fin-str)/CLOCKS_PER_SEC));

	str = clock();

        // Display

	if(iprint)
	  for(i = 0 ; i < nresIF ; i++)
	       printf("node %3d    %2d    %s \n", i , node_list[i] , RecLigAg->idres[res_list->ifres_num[i]] );

	// Applying the MWIS results into Atom Group Structure
        int res_i;
        int rot_i;

	#pragma omp parallel for private(res_i, rot_i)
        for( i = 0 ; i < nresIF ; i++){
                res_i = res_list->ifres_num[i];
                rot_i = node_list[i] ;
//              rot_i = node_list[i] ;
                apply_rot_to_atmgrp( RecLigAg , rotinf[i] , res_i , rot_i );
        }

	// Freeing
        Free_graph_t(&graph, nnodes, ncliques_mal);

	#pragma omp parallel for
	for(i = 0 ; i < nresIF ; i++)
		Free_rotinf(&rotinf[i]);

        fin = clock();
        if(iprint) printf("Part 4 : %f   sec\n", ((double)(fin-str)/CLOCKS_PER_SEC));


}

static void pack_reslist_1( struct atomgrp* RecLigAg , struct ifres_list* res_list , struct rot_info* rotinf_in, int* indexN)
{
	int i,j,k,l;

        int nresIF = res_list->num_of_ifres;

        struct node_energy* energy;
        struct link_list* ln_list;
        node_energy_malloc( &energy , nresIF );
        link_list_malloc( &ln_list, nresIF );

        struct rot_info* rotinf[nresIF];
	#pragma omp parallel for
	for(i = 0 ; i < nresIF ; i++)
                rot_info_malloc(&rotinf[i]);

        int i_total;

	#pragma omp parallel for private(k,l)
        for(i = 0 ; i < nresIF ; i++){
                i_total = indexN[i];
                rotinf[i]->res_num = rotinf_in[i_total].res_num;
                rotinf[i]->natoms  = rotinf_in[i_total].natoms;
                rotinf[i]->nrotamers = rotinf_in[i_total].nrotamers;
                for(k = 0 ; k <= rotinf[i]->nrotamers ; k++){
                        rotinf[i]->rot[k].probability = rotinf_in[i_total].rot[k].probability;
                        for(l = 0 ; l < rotinf[i]->natoms ; l++){
                                rotinf[i]->rot[k].X[l] = rotinf_in[i_total].rot[k].X[l];
                                rotinf[i]->rot[k].Y[l] = rotinf_in[i_total].rot[k].Y[l];
                                rotinf[i]->rot[k].Z[l] = rotinf_in[i_total].rot[k].Z[l];
                        }
                }
        }

        // Probability Correction

        double prob_sum;
	#pragma omp parallel for private(j)
        for(i = 0 ; i < res_list->num_of_ifres ; i++){

           if( rotinf[i]->nrotamers > 0 ){
               rotinf[i]->rot[0].probability = nb_rot_prob ;
               for(j = 1 ; j <= rotinf[i]->nrotamers ; j++)
                       prob_sum = prob_sum + rotinf[i]->rot[j].probability ;
               for(j = 1 ; j <= rotinf[i]->nrotamers ; j++)
                        rotinf[i]->rot[j].probability = rotinf[i]->rot[j].probability * (1 - nb_rot_prob ) / prob_sum ;
           }
           else {
                rotinf[i]->rot[0].probability = 1 ;
           }
        }
        nresIF = res_list->num_of_ifres;

        // Generating Energy Graph

        eng_linklist_gen( RecLigAg , res_list , energy , ln_list , rotinf );    // OK

        int resi = 0;   // for consistency of notation

        // Finding the heaviest node

        double min_energy = energy->residues[resi][resi].rotamers[0][0];
        int min_rot = 0;

        for(int ri = 1 ; ri < rotinf[resi]->nrotamers ; ri++)
                if( energy->residues[resi][resi].rotamers[ri][ri] <  min_energy ){
                        min_energy = energy->residues[resi][resi].rotamers[ri][ri];
                        min_rot = ri;
                }

        int res_i = res_list->ifres_num[resi];

        if(iprint)
	        printf("node %3d      %s \n", 0 , RecLigAg->idres[res_i] );

        apply_rot_to_atmgrp( RecLigAg , rotinf[resi] , res_i , min_rot );

	#pragma omp parallel for
        for(i = 0 ; i < nresIF ; i++)
                Free_rotinf(&rotinf[i]);

        Free_link_list ( &ln_list , nresIF );
        Free_node_energy( &energy , nresIF );
}

static void pack_reslist_2( struct atomgrp* RecLigAg , struct ifres_list* res_list , struct rot_info* rotinf_in, int* indexN)
{
	int i,j,k,l;
	int nnodes;
        int nresIF = res_list->num_of_ifres;

        struct node_energy* energy;
        struct link_list* ln_list;
        node_energy_malloc( &energy , nresIF );
        link_list_malloc( &ln_list, nresIF );

        struct rot_info* rotinf[nresIF];
        for(i = 0 ; i < nresIF ; i++)
                rot_info_malloc(&rotinf[i]);

        int i_total;
	#pragma omp parallel for private(k,l)
        for(i = 0 ; i < nresIF ; i++){
                i_total = indexN[i];
                rotinf[i]->res_num = rotinf_in[i_total].res_num;
                rotinf[i]->natoms  = rotinf_in[i_total].natoms;
		rotinf[i]->nrotamers = rotinf_in[i_total].nrotamers;	
                for(k = 0 ; k <= rotinf[i]->nrotamers ; k++){
                        rotinf[i]->rot[k].probability = rotinf_in[i_total].rot[k].probability;
                        for(l = 0 ; l < rotinf[i]->natoms ; l++){
                                rotinf[i]->rot[k].X[l] = rotinf_in[i_total].rot[k].X[l];
                                rotinf[i]->rot[k].Y[l] = rotinf_in[i_total].rot[k].Y[l];
                                rotinf[i]->rot[k].Z[l] = rotinf_in[i_total].rot[k].Z[l];
                        }
                }
        }

	// Probability Correction

        double prob_sum;

	#pragma omp parallel for private(j)

        for(i = 0 ; i < res_list->num_of_ifres ; i++){

           if( rotinf[i]->nrotamers > 0 ){
               rotinf[i]->rot[0].probability = nb_rot_prob ;
               for(j = 1 ; j <= rotinf[i]->nrotamers ; j++)
                       prob_sum = prob_sum + rotinf[i]->rot[j].probability ;
               for(j = 1 ; j <= rotinf[i]->nrotamers ; j++)
                        rotinf[i]->rot[j].probability = rotinf[i]->rot[j].probability * (1 - nb_rot_prob ) / prob_sum ;
           }
           else {
                rotinf[i]->rot[0].probability = 1 ;
           }
        }

	reslist_nnodes( RecLigAg , res_list , &nnodes , energy , ln_list , rotinf );

	// Graph Generation Based on Res_List

        struct graph_t* graph ;
        graph_t_malloc( &graph , nnodes , 0);
        int node_info[nnodes][4];               // DICTIONARY : for node numbers conversion into rotamer data
        graph_construction( graph , res_list , ln_list , energy , rotinf , node_info , nnodes); 
        weight_edit(graph);             // Scale MAX - W(i) & Normalize W(i)'s to 1     
        Free_link_list ( &ln_list , nresIF );
        Free_node_energy( &energy , nresIF );

	// Finding the heaviest triple of nodes by search
	
	double max_weight = 0;
	int max_node = 0;

	#pragma omp parallel for
	for(i = 0 ; i < graph->vertex_count ; i++)
		if(graph->weight[i] > max_weight){
			max_weight = graph->weight[i];
			max_node = i;
		}

	// Finding the indices of the found nodes
	
	int node_list[nresIF];

	node_list[node_info[max_node][0]] = node_info[max_node][2]  ;
	node_list[node_info[max_node][1]] = node_info[max_node][3]  ;

        // Display

	if(iprint)
           for(i = 0 ; i < nresIF ; i++)
                printf("node %3d    %2d    %s \n", i , node_list[i] , RecLigAg->idres[res_list->ifres_num[i]] );

	// Applying the MWIS results into Atom Group Structure
	
        int res_i;
        int rot_i;

        for( i = 0 ; i < nresIF ; i++){
                res_i = res_list->ifres_num[i];
                rot_i = node_list[i] ;
                rot_i = node_list[i] ;
                apply_rot_to_atmgrp( RecLigAg , rotinf[i] , res_i , rot_i );
        }
  
	// Freeing

	#pragma omp parallel for
        for(i = 0 ; i < nresIF ; i++)
                Free_rotinf(&rotinf[i]);

        Free_graph_t(&graph, nnodes, 0);
}

static void index_fixer( struct ifres_list* res_list , struct ifres_list* resList_C[] , int nclusters , int** index_C )
{
        for(int i = 0 ; i < nclusters ; i++)
               for(int j = 0 ; j < resList_C[i]->num_of_ifres ; j++)
                       for( int r = 0 ; r < res_list->num_of_ifres ; r++){
                               if( resList_C[i]->ifres_num[j] == res_list->ifres_num[r]){
                                        index_C[i][j] = r ;
                                        break;
                                }
                        }
}

void full_pack( struct atomgrp* RecLigAg, struct List lig_list, struct rot_info* Major_rotinf, int nIFres , struct ifres_list* res_list)
{
	
	// Ligand-Receptor Separation

	struct List List_Rec;
	List_Rec.K = malloc( nIFres * sizeof(int));
	struct List List_Lig;
	List_Lig.K = malloc( nIFres * sizeof(int));
	int ind_r = 0 ;
        int ind_l = 0 ;

	int i,k,l;

	for(i = 0 ; i < nIFres ; i++){
		if( RecLigAg->iares[Major_rotinf[i].res_num] < lig_list.K[0] ){
			List_Rec.K[ind_r] = i ;
			ind_r++;
		}
		else{
			List_Lig.K[ind_l] = i ;
			ind_l++;
		}
	}
	List_Rec.n = ind_r ;
	List_Lig.n = ind_l ;

        int cluster_i , nclusters;

        struct rot_info* Minor_rotinf;	
        mal_rotinf(nIFres , &Minor_rotinf);	

        IFres_adjustment( RecLigAg, res_list, &cluster_i, Major_rotinf, nIFres, Minor_rotinf );

	// continue...

	if( res_list->num_of_ifres != 0 ){

		struct ifres_list* resList_C[res_list->num_of_ifres];
                for(i = 0 ; i < res_list->num_of_ifres ; i++)
                        ifres_Clus_malloc( &resList_C[i], res_list->num_of_ifres);

                reslist_clustering( RecLigAg , res_list , resList_C , &nclusters);

		int nresIF = res_list->num_of_ifres;

	 	if(iprint) printf("Nres = %d     Nclus = %d\n",nresIF, nclusters);

		if(UNB_i)	// To Exclude The Unbound Conformers	
		        for(i = 0 ; i < nresIF ; i++){
                           if(Minor_rotinf[i].nrotamers > 1){  				   
				Minor_rotinf[i].nrotamers = Minor_rotinf[i].nrotamers - 1;
		                for(k = 0 ; k <= Minor_rotinf[i].nrotamers ; k++){
		                        Minor_rotinf[i].rot[k].probability = Minor_rotinf[i].rot[k+1].probability;
		                        for(l = 0 ; l < Minor_rotinf[i].natoms ; l++){
		                                Minor_rotinf[i].rot[k].X[l] = Minor_rotinf[i].rot[k+1].X[l];
                                                Minor_rotinf[i].rot[k].Y[l] = Minor_rotinf[i].rot[k+1].Y[l];
                                                Minor_rotinf[i].rot[k].Z[l] = Minor_rotinf[i].rot[k+1].Z[l];
					}
				}
			    }
			}	


                int** index_C;
                index_C = malloc ( nclusters * sizeof(int*));
                for(i = 0 ; i < nclusters ; i++)
                        index_C[i] = malloc ( resList_C[i]->num_of_ifres * sizeof(int));

                index_fixer( res_list , resList_C , nclusters , index_C );


	        //clock_t start, finish;
	        //start = clock();

		clock_t str, fin;
	
		#pragma omp parallel for //private(RecLigAg)	// def : not commented out
                for(i = 0 ; i < nclusters ; i++){
			str = clock();						

                        if(resList_C[i]->num_of_ifres == 1)
                                pack_reslist_1( RecLigAg , resList_C[i] , Minor_rotinf , index_C[i] );
			else if(resList_C[i]->num_of_ifres == 2)
				pack_reslist_2( RecLigAg , resList_C[i] , Minor_rotinf , index_C[i] );
			else
	                        pack_reslist( RecLigAg , resList_C[i] , Minor_rotinf , index_C[i] );

		        fin = clock();
			if( (0) && (resList_C[i]->num_of_ifres > 2))  
				printf("Cluster = %d  with nRes = %d  :  Packing Time = %f\n",i, resList_C[i]->num_of_ifres, 							((double)(fin-str)/CLOCKS_PER_SEC));
		}

        	//finish = clock();
		//	if(1) printf("Packing Time = %f\n",((double)(finish-start)/CLOCKS_PER_SEC));

	        //	if(1) printf("%d      %d	%f\n", nresIF, nclusters, ((double)(finish-start)/CLOCKS_PER_SEC));

                for(i = 0 ; i < res_list->num_of_ifres ; i++)
                        Free_ifres_list(&resList_C[i]);
                for(i = 0 ; i < nclusters ; i++)
                        free(index_C[i]);
                free(index_C);
		free(List_Rec.K);
		free(List_Lig.K);
	        Free_rotinf(&Minor_rotinf);
    }
}
