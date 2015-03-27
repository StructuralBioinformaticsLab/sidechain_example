
#ifndef _ROTAMER_H
#define _ROTAMER_H

/*  INTER_RES_H_ */

struct ifres_list
{
        int num_of_ifres;
        int* ifres_num ;
};

/* _ROTAMER_PLACEMENT_H */

struct dunbrack_rotamer {

        int rotamer_number;
        double* probability;
        int16_t* chi1;
        int16_t* chi2;
        int16_t* chi3;
        int16_t* chi4;
};

void ifres_list_malloc( struct ifres_list** res_list );
void ifres_Clus_malloc( struct ifres_list** res_list , int nres_IF );
void Free_ifres_list( struct ifres_list** res_list );
void array_sort( int* array , int size );
void IFres_adjustment(struct atomgrp* RecLigAg, struct ifres_list* res_list, int* cluster_i, struct rot_info* Major_rotinf, int nIFres , struct rot_info* Minor_rotinf);
void rot_info_malloc(struct rot_info** rotinf);
void Free_rotinf(struct rot_info** rotinf);
int load_file_to_memory(const char *lib_filename, char **lib_content);
void get_rotamer_coordinates( struct atomgrp* ag,  int resi_num, char* lib_contents, struct rot_info* rotinf, int number_rotamers  );
void apply_rot_to_atmgrp( struct atomgrp* ag,  struct rot_info* rotinf, int resi_number, int rotamer_index );

#endif
