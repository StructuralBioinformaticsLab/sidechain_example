#ifndef _GRAPH_H
#define _GRAPH_H

struct graph_t {
	long vertex_count;
	double *weight;
	unsigned char *adjacency_matrix;

	int ncliques;	
	int** clique_list;	// clique_list[i] : list of nodes of clique i
	int*  clique_size;	// clique_size[i] : size of clique i
	int** cliquesof_i;	// cliquesof_i[i] : list of cliques (and edges) which include node i
	int*  node_nClq;	// node_nClq[i]   : number of cliques which include node i
};

//      link_list->residues[i][j].rotamers[r][k] = 0/1  => 0 : rotamer i_r and j_k are dead-ending and eliminated from active list
struct link_list {
        int  nresidue ;
        int* nrotamers;
        struct rotamer_pair** residues;
};

struct rotamer_pair{
        int** rotamers; // 0-1 indicates whether the pair i_r & j_k is linked or not
};

//      node_energy->residues[i][j].rotamers[r][k] => the value of the interaction energy between rotamer i_r and j_k
struct node_energy{
        int  nresidue ;
        int* nrotamers;
        struct eng_residue** residues;
};

struct eng_residue{
        double** rotamers; // Inter_Energy of rotamer pair
};



#define graph_edge(graph, i, j)    ((graph)->adjacency_matrix[(i) * (graph)->vertex_count + (j)])

void graph_t_malloc( struct graph_t** graph , int MAX_NODE, int N_CLIQUES);
void Free_graph_t  ( struct graph_t** graph , int MAX_NODE, int N_CLIQUES);
void weight_edit(struct graph_t* graph);
void eng_linklist_gen( struct atomgrp *RecLigAg , struct ifres_list* res_list , struct node_energy* energy , struct link_list* ln_list ,                          struct rot_info* rotinf[]);
void node_counter( struct ifres_list* res_list , struct link_list* ln_list , struct rot_info* rotinf[], int* nnodes);
void graph_construction( struct graph_t* graph , struct ifres_list* res_list , struct link_list* ln_list ,  struct node_energy* energy , 			    struct rot_info* rotinf[] , int node_info[][4] , int nnodes);
void add_2cliques_to_graph( struct graph_t* graph);


//#ifndef _DEE_H

void dee( struct node_energy* energy , struct link_list* ln_list );
void link_list_malloc( struct link_list** ln_list , int MAX_RES1 );
void Free_link_list ( struct link_list** ln_list , int MAX_RES1 );
void node_energy_malloc( struct node_energy** energy , int MAX_RES1 );
void Free_node_energy( struct node_energy** energy , int MAX_RES1 );

void vdweng_inter(struct atomgrp *ag, double rc, int *list1, int numOfAtomsList1, int *list2, int numOfAtomsList2, double *eng);
void get_list( struct atomgrp *ag , int resid , int* list, int* list_size);

#endif
