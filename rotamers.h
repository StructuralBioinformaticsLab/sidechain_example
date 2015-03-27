#ifndef ROTAMERS_H
#define ROTAMERS_H

#define MAX_NUM_OF_GROUPS_INTERFACE   100
#define MAX_NUM_OF_RES_INTERFACE      200
#define ENERGY_THRESHOLD_REPACK       10// kcal/mol
#define PROBABILITY_THRESHOLD         0.99

/* Defining the Structs to store Rotamers' Information & Coordinates */

struct res_groups {
  int res_num[200];
  int num_of_res;
};

struct adjacency_list {
  int node;
  int next;
};

void init_rotinf(struct atomgrp *ag, int num_of_res_interface, int *res_list_interface, char *lib_contents, struct rot_info **rotinf);

void mal_rotinf(int num_of_res_interface, struct rot_info **rotinf);

#endif
