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
#include <omp.h>
#include _MOL_INCLUDE_
#include "main.h"
#include "parameters.h"
#include "rotamer_placement.h"
#include "pack.h"
#include "rotamers.h"
#include "utils.h"


int main(int argc, char* argv[]){

	// Input Files //
	/* Intro Message */
                printf("\n");
                printf("A Program to Repack Protein Side-chains for Protein Docking Refinement Procedures\n");
                printf("Copyright (c) 2014, Structural Bioinformatics Laboratory, Boston University\n");
                printf("Author: Mohammad Moghadasi (mohamad@bu.edu) \n");	

	if(argc!=10){
		printf("Usage:\n ./main \n   Complex_IN.pdb Complex_IN.psf Complex_IN.mol2 Complex_IN_Ligand.pdb\n   Libmol-param-file charmm-param-file.prm  charmm-rtf-file.rtf rotamer-library-binary-file.txt\n   Complex_OUT.pdb\n \n");
		exit(EXIT_FAILURE);
	}

        char* ifile           = argv[1];       //pdb file of both receptor and ligand
        char* psffile         = argv[2];         //charmm type psf file
        char* mol2file        = argv[3];         //mol2 file
        char* pdbfilelig      = argv[4];         //pdb file of ligand
        char* atom_prm_file   = argv[5];         //libmol parameter file
              prmfile         = argv[6];         //charmm type parameter file
              rtffile         = argv[7];         //charmm type connectivity file
        char* rotamer_library_file  = argv[8]; //rotamer raw library file
        char* ofile           = argv[9];         //output file
	
	// Filling the atom_group struct //

	struct prm *atomprm    = read_prm(atom_prm_file,_MOL_VERSION_);
	struct atomgrp* ag     = read_file_atomgrp(ifile, atomprm, -1);
	read_ff_charmm(psffile, prmfile, rtffile, ag);

	if(!read_hybridization_states_from_mol2(mol2file,ag)){
	    exit (EXIT_FAILURE);             
	}                 
        fix_acceptor_bases(ag,atomprm);

        struct List lig_list;
	read_fix(pdbfilelig,&lig_list.n,&lig_list.K);

	fixed_init(ag);
        fixed_update_unfreeze_all(ag);
        zero_grads(ag);
        fill_ingrp(ag);

	struct agsetup* ags;
        ags     = malloc(sizeof(struct agsetup));
        init_nblst(ag,ags);
        update_nblst(ag,ags);

	// Mark interface residues //
	
        int num_of_res_interface;
        int res_list_interface[ag->nres];
        mark_interface_residues(ag,ags,lig_list, lig_rec_dist ,&num_of_res_interface,res_list_interface);

	// Initialize side chain rotamer library  //
	
	//nrotCoef = 3;
	nrotCoef = 1;
	MAX_ROT = 245;
	MAX_RES = ag->nres;//needed for full_pack
        cutoff = 3;

	// Reslist // 

	struct ifres_list* reslist;
	ifres_list_malloc( &reslist );
	reslist->num_of_ifres = num_of_res_interface;
	for(int r = 0; r < reslist->num_of_ifres ; r++)
		reslist->ifres_num[r] = res_list_interface[r];

	// Library //

	char *rotamer_lib;
	load_file_to_memory(rotamer_library_file, &rotamer_lib);
	struct rot_info *rotinf;
	init_rotinf(ag, num_of_res_interface, res_list_interface, rotamer_lib, &rotinf);

	struct ifres_list* reslist_minor;
	ifres_list_malloc( &reslist_minor ) ;

	// MAIN FUNCITION // 
	//
	clock_t start, finish;
	start = clock();

	full_pack(ag,lig_list,rotinf,num_of_res_interface, reslist_minor);

	finish = clock();
	if(0) printf("Processing Time = %f\n",((double)(finish-start)/CLOCKS_PER_SEC));

	// Writing the atom_group into a PDB //

	write_pdb_traj_nopar(ag,ifile,ofile);
	
	// Free memory //
        Free_ifres_list( &reslist_minor );	
	free(rotamer_lib);

	return 0;
}
