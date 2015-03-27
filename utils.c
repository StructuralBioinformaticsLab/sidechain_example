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
#include _MOL_INCLUDE_
#include "main.h"
#include "utils.h"


void read_fix(char *ffile, int *nfix, int **fix){
   int linesz=91;
   char *buffer=malloc(sizeof(char)*linesz);
   *nfix=0;
   FILE* fp = fopen (ffile, "r");
   while (fgets(buffer, linesz-1, fp)!=NULL)
   {
       if(!strncmp(buffer,"ATOM",4))(*nfix)++;
   }
   fclose(fp);
   *fix=malloc(*nfix*sizeof(int));
   fp = fopen (ffile, "r");
   int na=0;
   while(fgets(buffer, linesz-1, fp)!=NULL)
   {
       if(!strncmp(buffer,"ATOM",4))
       {
         (*fix)[na]=atoi(buffer+4)-1;
         na++;
       }
   }
   free(buffer);
 //  free(fix);
   fclose(fp);
}


void fixed_update_unfreeze_all(struct atomgrp *ag){

	int i, m, n;
	struct atom *a0, *a1, *a2, *a3;
	if (ag->nbact > 0) {
		free(ag->bact);
		ag->nbact = 0;
	}
	if (ag->nangact > 0) {
		free(ag->angact);
		ag->nangact = 0;
	}

	if (ag->ntoract > 0) {
		free(ag->toract);
		ag->ntoract = 0;
	}
	if (ag->nimpact > 0) {
		free(ag->impact);
		ag->nimpact = 0;
	}
	if (ag->nactives > 0) {
		free(ag->activelist);
		ag->nactives = 0;
	}
// atoms
	for (i = 0; i < ag->natoms; i++)
		ag->atoms[i].fixed = 0;

	int ci;
	ag->activelist = (int *)_mol_malloc((ag->natoms) * sizeof(int));
	ci = 0;
	for (i = 0; i < ag->natoms; i++) {
		if (ag->atoms[i].fixed == 0)
			ag->activelist[ci++] = i;
	}
	ag->nactives = ci;
// bonds
	n = ag->nbonds;
	m = 0;
	ag->bact = _mol_malloc(n * sizeof(struct atombond *));
	for (i = 0; i < n; i++) {
		a0 = ag->bonds[i].a0;
		a1 = ag->bonds[i].a1;
		if (a0->fixed == 1 && a1->fixed == 1)
			continue;
		ag->bact[m++] = &(ag->bonds[i]);
	}
	ag->nbact = m;
	if (m > 0)
		ag->bact =
		    (struct atombond **)_mol_realloc(ag->bact,
						     m *
						     sizeof(struct atombond *));
	else
		free(ag->bact);
// angles
	n = ag->nangs;
	m = 0;
	ag->angact = _mol_malloc(n * sizeof(struct atomangle *));
	for (i = 0; i < n; i++) {
		a0 = ag->angs[i].a0;
		a1 = ag->angs[i].a1;
		a2 = ag->angs[i].a2;
		if (a0->fixed == 1 && a1->fixed == 1 && a2->fixed == 1)
			continue;
		ag->angact[m++] = &(ag->angs[i]);
	}
	ag->nangact = m;
	if (m > 0)
		ag->angact =
		    (struct atomangle **)_mol_realloc(ag->angact,
						      m *
						      sizeof(struct atomangle
							     *));
	else
		free(ag->angact);
// dihedrals
	n = ag->ntors;
	m = 0;
	ag->toract = _mol_malloc(n * sizeof(struct atomtorsion *));
	for (i = 0; i < n; i++) {
		a0 = ag->tors[i].a0;
		a1 = ag->tors[i].a1;
		a2 = ag->tors[i].a2;
		a3 = ag->tors[i].a3;
		if (a0->fixed == 1 && a1->fixed == 1 && a2->fixed == 1
		    && a3->fixed == 1)
			continue;
		ag->toract[m++] = &(ag->tors[i]);
	}
	ag->ntoract = m;
	if (m > 0)
		ag->toract =
		    (struct atomtorsion **)_mol_realloc(ag->toract,
							m *
							sizeof(struct
							       atomtorsion *));
	else
		free(ag->toract);
// impropers
	n = ag->nimps;
	m = 0;
	ag->impact = _mol_malloc(n * sizeof(struct atomimproper *));
	for (i = 0; i < n; i++) {
		a0 = ag->imps[i].a0;
		a1 = ag->imps[i].a1;
		a2 = ag->imps[i].a2;
		a3 = ag->imps[i].a3;
		if (a0->fixed == 1 && a1->fixed == 1 && a2->fixed == 1
		    && a3->fixed == 1)
			continue;
		ag->impact[m++] = &(ag->imps[i]);
	}
	ag->nimpact = m;
	if (m > 0)
		ag->impact =
		    (struct atomimproper **)_mol_realloc(ag->impact,
							 m *
							 sizeof(struct
								atomimproper
								*));
	else
		free(ag->impact);
}


void mark_interface_residues(struct atomgrp* ag, struct agsetup* ags, struct List lig_list, double cutoff, int *num_of_interface_res, int *interface_res){
  int i,j;
  int ii,jj;
  int flag;
  int first_atom, end_atom;
  int first_ligand_atom;
  int markag[ag->natoms];
  double dx,dy,dz,r2;
  double cutoffSq = cutoff*cutoff; 

  first_ligand_atom = lig_list.K[0];
  //printf("first ligand atom: %d\n",first_ligand_atom);
  //printf("nfat: %d\n",ags->nblst->nfat);
  
  //set markag array to zero
  for(i=0;i<ag->natoms;i++){
    markag[i] = 0;    
  }

  //mark atoms within cutoff with 1
  for(i=0;i<ags->nblst->nfat;i++){
    ii  = ags->nblst->ifat[i];
    //printf("ii: %d\n",ii);

    for(j=0; j<ags->nblst->nsat[i]; j++){
      jj = ags->nblst->isat[i][j];
      //printf("ligand ii: %d, receptor jj: %d\n",ii,jj);
	
      if((ii<first_ligand_atom && jj>=first_ligand_atom) || (ii>=first_ligand_atom && jj<first_ligand_atom)){//ligand atoms
	  //printf("ligand ii: %d, receptor jj: %d\n",ii,jj);
	  
	dx = ag->atoms[ii].X - ag->atoms[jj].X;
	dy = ag->atoms[ii].Y - ag->atoms[jj].Y;
	dz = ag->atoms[ii].Z - ag->atoms[jj].Z;
	r2 = dx*dx + dy*dy + dz*dz;
    
	if(r2 <= cutoffSq){
	  markag[ii] = 1;
	  markag[jj] = 1;
	  //printf("atoms %d and %d are withing cutoff range\n",ii,jj);
	}
      }
    }
  }

  ii = 0;
  for(i=0;i<ag->nres;i++){
    first_atom = ag->iares[i];    
    if(i<ag->nres-1){ 
      end_atom = ag->iares[i+1]; 
    }else{
      end_atom = ag->natoms;
    }
    
    flag = 0;
    for(j=first_atom;j<end_atom;j++){
      if(markag[j]==1){ 
	flag = 1;
	break;
      }
    }

    if(flag==1){
      interface_res[ii]=i;
      ii++;
      //printf("   test residue: %d\n",i);
    }
  }

  (*num_of_interface_res) = ii;
}


void write_pdb_traj_nopar(struct atomgrp* ag, const char* inf, const char* ouf){
        FILE* fp = myfopen (inf, "r");
        FILE* fop = myfopen(ouf, "w");

        char* line = NULL;
        size_t len = 0;
        char c;

        //write MODEL
        //fprintf(fop,"MODEL\n");

        // read every line of the pdb file
        int atomi = 0;
        while (getline (&line, &len, fp) != -1)
        {
                if (strncmp (line, "ATOM  ", 6) != 0 && strncmp (line, "HETATM", 6) != 0)
                {
                         fprintf (fop,"%s",line);
                        continue;
                }
                c=line[54];
                sprintf(line+30,"%8.3f",ag->atoms[atomi].X );
                sprintf(line+38,"%8.3f",ag->atoms[atomi].Y );
                sprintf(line+46,"%8.3f",ag->atoms[atomi].Z );
                atomi++;
                line[54]=c;
                fprintf (fop,"%s",line);
        }
        //fprintf(fop,"ENDMDL\n");

        if (line)
                free (line);
        myfclose (fp);
        myfclose (fop);
}

