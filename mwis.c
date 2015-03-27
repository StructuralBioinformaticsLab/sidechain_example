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
#include <limits.h>
#include <stdint.h>
#include <time.h>
#include _MOL_INCLUDE_
#include "main.h"
#include "parameters.h"
#include "rotamer_placement.h"
#include "graph.h"
#include "omp.h"
#include "mwis.h"

#define min(x, y)    ((x) < (y) ? (x) : (y))
#define max(x, y)    ((x) > (y) ? (x) : (y))

static void y_i_theta_singh (struct graph_t *graph , double* theta, int i, double epsilon ,double* y_i_theta)
{

        // Calculating a_i(theta) 
        double sum_theta = 0;

        for(int j = 0 ; j < graph->node_nClq[i]; j++)	// all cliques of node i 
                sum_theta += theta[j];


        // Calculating x_i(theta)

        double y_i_t;

        if(sum_theta >= graph->weight[i])
                y_i_t = 0;
        else if (sum_theta <= graph->weight[i] - epsilon)
                y_i_t = 1;
        else
                y_i_t = ( graph->weight[i] - sum_theta ) / epsilon;

        y_i_theta[0] = y_i_t;
}

static void x_i_theta_calc (struct graph_t *graph , double * theta, int i, double epsilon ,double* x_i_theta)
{
	// Calculating a_i(theta) 
	double a_i_theta = graph->weight[i];

	for(int j = 0 ; j < graph->node_nClq[i]; j++)
		a_i_theta = a_i_theta - theta[graph->cliquesof_i[i][j]];

	// Calculating x_i(theta)
	
	double x_i_t;

	if(a_i_theta > 0)
		x_i_t = ( 1 - 2*epsilon/a_i_theta + sqrt(4*pow(epsilon/a_i_theta,2) + 1) ) / 2;
	else if (a_i_theta < 0)
		x_i_t = ( 1 - 2*epsilon/a_i_theta - sqrt(4*pow(epsilon/a_i_theta,2) + 1) ) / 2;
	else 
		x_i_t = 0.5;

	x_i_theta[0] = x_i_t;
}

void gradproj_CC (struct graph_t *graph , double * x_gp)
{
	long i,j,k,n;

        clock_t str, fin;

	// Step 1. Initialization
	
	// I.

	int N = N_mw1;
        double lambda = lambda_mw1;
        double epsilon = epsilon_mw1;
        double eps_lowbound = eps_lowbound_mw1;

	double theta_new[graph->ncliques];
	double theta_old[graph->ncliques];
	double x_new[graph->vertex_count];
	double x_old[graph->vertex_count];

        str = clock();

        #pragma omp parallel for private(i)
	for(j = 0 ; j < graph->ncliques ; j++){
		theta_old[j] = graph->weight[graph->clique_list[j][0]];
		for(i = 0 ; i < graph->clique_size[j] ; i++)	// all nodes in clique j
			if(graph->weight[graph->clique_list[j][i]] > theta_old[j])
				theta_old[j]= graph->weight[graph->clique_list[j][i]];
	}

        fin = clock();
        if(0) printf("MWIS 1  : %f   sec\n", ((double)(fin-str)/CLOCKS_PER_SEC));

	str = clock();


	// II.
	double x_i_theta;
        
	#pragma omp parallel for private(x_i_theta)
	for(i = 0 ; i < graph->vertex_count ; i++){
		x_i_theta_calc (graph , theta_old, i, epsilon, &x_i_theta);
		x_old[i] = x_i_theta;
	}
	
	fin = clock();
        if(0) printf("MWIS 2  : %f   sec\n", ((double)(fin-str)/CLOCKS_PER_SEC));

	str = clock();


	// Step 2.	
//	int stop_loop = 0; 	// 0 : continue
	n = 1;			// n : iteration index

//	while ( stop_loop == 0 || n < N ){

//	#pragma omp parallel for	
	for(n = 1; n < N; n++){
//	while(n<N){

		// Part b

		epsilon = (epsilon < eps_lowbound) ? eps_lowbound : epsilon * 0.9;
		 
   	        #pragma omp parallel for private(k)
		for(j = 0 ; j < graph->ncliques ; j++){
			double sum = 0;
//			#pragma omp parallel for ordered 
			for(k = 0 ; k < graph->clique_size[j] ; k++) 
				sum = sum + x_old[graph->clique_list[j][k]];
			theta_new[j] = theta_old[j] - lambda*(1-sum);			
			if( theta_new[j] < 0 )
				theta_new[j] = 0 ;
		}	
	    
                #pragma omp parallel for private(x_i_theta)
		for(i = 0 ; i < graph->vertex_count ; i++){
			x_i_theta_calc (graph , theta_new, i, epsilon, &x_i_theta);
			x_new[i] = x_i_theta;
		}
	    
		// Stopping Criteria
		/*
		stop_loop = 1;

		for(j = 0 ; j < graph->ncliques ; j++)
			if(theta_new[j] - theta_old[j] > delta_mw1){
				stop_loop = 0;
				break;
			}	
		*/
			
		//	n++;

                #pragma omp parallel for		
		for(j = 0 ; j < graph->ncliques ; j++)	theta_old[j] = theta_new[j];
	        #pragma omp parallel for	
		for(i = 0 ; i < graph->vertex_count ; i++)	x_old[i] = x_new[i];
	}

	fin = clock();
        if(0) printf("MWIS 3  : %f   sec\n", ((double)(fin-str)/CLOCKS_PER_SEC));

	str = clock();

	#pragma omp parallel for
	for(i = 0 ; i < graph->vertex_count ; i++)
		x_gp[i] = x_new[i];

	fin = clock();
        if(0) printf("MWIS 4  : %f   sec\n", ((double)(fin-str)/CLOCKS_PER_SEC));
}


void myestimation(struct graph_t *graph, double *lp_x, unsigned char *x)
{
        unsigned char *x1;
        long i, j;
        const unsigned long n = graph->vertex_count;

        if (!x)
                return;
        x1 = (unsigned char *) malloc(n);
        if (!x1)
        {
                return;
        }

        /* initialization */
        #pragma omp parallel for
        for (i = 0; i < graph->vertex_count; i++)
        {
                if (lp_x[i] >= 1 -  myTHRESHOLD1)      x[i] = 1;
                else if (lp_x[i] <= myTHRESHOLD0)      x[i] = 0;
                else                               x[i] = 2; /* undetermined */
        }

        // Set x[i] to 0 if i has a neighbor of value 1
	for (i = 0; i < graph->vertex_count; i++)
		for (j = 0; j < graph->vertex_count; j++)
			if( graph_edge(graph, i, j) && (x[j]==1)){
					x[i]=0;
					break;
				}

	// Set x[i] to 1 if all it's neighbors are 0
	int ind_1;
	for (i = 0; i < graph->vertex_count; i++){
		ind_1 = 1;
		for (j = 0; j < graph->vertex_count; j++){
			if( graph_edge(graph, i, j) && (x[j]!=0)){
				ind_1 = 0;
				break;
			}
		}				
		if(ind_1)
			x[i] = 1;	
	}

        // Undetermined values
	for (i = 0; i < graph->vertex_count; i++)
		if(x[i] == 2){
			x[i] = 1;
			for (j = 0; j < graph->vertex_count; j++)
				if( graph_edge(graph, i, j))
					x[j] = 0;
		}	
}

void gradproj_singh (struct graph_t *graph , double * x_gp)
{

	long i,j,k;

	int MAX_nCliq_i = -1;
	for(long i = 0 ; i < graph->vertex_count ; i++)
		MAX_nCliq_i = (MAX_nCliq_i < graph->node_nClq[i]) ? graph->node_nClq[i] : MAX_nCliq_i;

	double beta_old[graph->vertex_count];
	double beta_new[graph->vertex_count];
	double theta_old[graph->vertex_count][MAX_nCliq_i];
	double theta_new[graph->vertex_count][MAX_nCliq_i];
        double eta_new[graph->vertex_count][MAX_nCliq_i];		

	#pragma omp parallel for private(j)	
	for(i = 0 ; i < graph->vertex_count ; i++){
		beta_old[i] = 0;
                beta_new[i] = 0;
		#pragma omp parallel for
		for(j = 0 ; j < MAX_nCliq_i ; j++){
			theta_old[i][j] = 0;
                        theta_new[i][j] = 0;
                        eta_new[i][j]   = 0;

		}
	}	

	// Initialize

       #pragma omp parallel for private(j,k)
       for(i = 0 ; i < graph->vertex_count ; i++){
                beta_old[i] = 1;

		#pragma omp parallel for private(k)
                for(j = 0 ; j < graph->node_nClq[i] ; j++){
                        theta_old[i][j] = graph->weight[graph->clique_list[graph->cliquesof_i[i][j]][0]];

                        // all nodes in jth clique of node i			
			#pragma omp parallel for
                        for(k = 0 ; k < graph->clique_size[graph->cliquesof_i[i][j]] ; k++)
                                if(graph->weight[graph->clique_list[graph->cliquesof_i[i][j]][k]] > theta_old[i][j])
                                        theta_old[i][j]= graph->weight[graph->clique_list[graph->cliquesof_i[i][j]][k]];

                }

                // display
                if(0){
                        printf("graph->ncliques = %d\n",  graph->ncliques);
                        printf("graph->node_nClq[%ld] = %d\n", i,graph->node_nClq[i]);
                        for(j = 0 ; j < graph->node_nClq[i] ; j++){
                                printf("theta_old[%ld][%ld] = %f\n", i,j,theta_old[i][j]);
                                printf("graph->cliquesof_i[%ld][%ld] = %d\n", i,j,graph->cliquesof_i[i][j]);
                        }
                        printf("\n");
                }
	}

        // Repeat    
          
	int n;                
//	int stop_loop = 0;                

	int N = N_mw2;
	double epsilon = epsilon_mw2;
	double eps_lower_bound = eps_lowbound_mw2;
	
	//  	while ( stop_loop == 0 && n < N ){

	#pragma omp parallel for
        for(n = 1; n < N; n++){

		// n = n+1;

		epsilon = ( (epsilon / 2) < eps_lower_bound ) ? eps_lower_bound : epsilon / 2;

                double y_i_theta;
                double sum_y_j_eta;
                double sum_c_j;

		#pragma omp parallel for private(j,k, sum_c_j, sum_y_j_eta, y_i_theta)
	        for(i = 0 ; i < graph->vertex_count ; i++){
        
                        beta_new[i] = ( 1 + pow((1+4*pow(beta_old[i],2)),0.5)) / 2;
                        
                        if(0)
                           printf("beta_new[%ld] = %f\n",i,beta_new[i]);

                        if(0)
                           printf("For size: graph->node_nClq[%ld] = %d\n",i,graph->node_nClq[i]);

			int kk;                        
			#pragma omp parallel for private(sum_c_j, sum_y_j_eta, y_i_theta, k, kk)
                        for(j = 0 ; j < graph->node_nClq[i] ; j++){         // all cliques of node i
                        
                                sum_y_j_eta = 0;
                                sum_c_j = 0;

                                for(k = 0 ; k < graph->clique_size[graph->cliquesof_i[i][j]] ; k++){// all nodes of jth clq of ith node

                                        kk = graph->clique_list[graph->cliquesof_i[i][j]][k];       // actual node number
        
                                        y_i_theta_singh (graph , eta_new[kk], kk, epsilon , &y_i_theta);

                                        sum_y_j_eta += y_i_theta;

                                        sum_c_j += graph->node_nClq[kk];
                                }

                                theta_new[i][j]=  eta_new[i][j] - epsilon * (1 - sum_y_j_eta) / sum_c_j ; 
                                theta_new[i][j] = ( theta_new[i][j] > 0 ) ? theta_new[i][j] : 0;
                        }
                        
			#pragma omp parallel for
                        for(j = 0 ; j < graph->node_nClq[i] ; j++){         // all cliques of node i

                                eta_new[i][j] = theta_new[i][j] + ((beta_old[i] - 1)/beta_new[i])*(theta_new[i][j] - theta_old[i][j]);
                                
                        }

			#pragma omp parallel for
                        for(j = 0 ; j < graph->node_nClq[i] ; j++){
                                theta_old[i][j] = theta_new[i][j];
                        }     
   		} 
    	}     

	// Extract Primal Variables Values
	
	double y_i_theta;	

	#pragma omp parallel for private(y_i_theta)
	for(i = 0 ; i < graph->vertex_count ; i++){
		y_i_theta_singh (graph , theta_old[i], i, epsilon , &y_i_theta);
		x_gp[i] = y_i_theta;
	}

	// Display
	if(0)
	   for(i = 0 ; i < graph->vertex_count ; i++)
		printf("Node %ld : weight = %f \n",i,graph->weight[i]);
}
