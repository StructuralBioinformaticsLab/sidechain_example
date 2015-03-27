#ifndef UTILS_H
#define UTILS_H

void read_fix(char *ffile, int *nfix, int **fix);

void fixed_update_unfreeze_all(struct atomgrp *ag);

void mark_interface_residues(struct atomgrp* ag, struct agsetup* ags, struct List lig_list, double cutoff, int *num_of_interface_res, int *interface_res);

void write_pdb_traj_nopar(struct atomgrp* ag, const char* inf, const char* ouf);

#endif
