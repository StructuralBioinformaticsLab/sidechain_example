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


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include _MOL_INCLUDE_
#include "main.h"
#include "parameters.h"
#include "utils.h"
#include "rotamers.h"
#include "rotamer_placement.h"
#include "utils.h"

void mal_rotinf(int num_of_res_interface, struct rot_info **rotinf){
  int i,j;

  (*rotinf) = (struct rot_info *) _mol_malloc(num_of_res_interface*sizeof(struct rot_info));

  for(i=0;i<num_of_res_interface;i++){

    (*rotinf)[i].rot = (struct rot_list *)_mol_malloc(MAX_ROT * sizeof(struct rot_list));

    float* memX;
    float* memY;
    float* memZ;

    for(j = 0 ; j < MAX_ROT ; j++){

      memX =  _mol_malloc (MAX_ATOM * sizeof(float));
      (*rotinf)[i].rot[j].X = memX;
      memY =  _mol_malloc (MAX_ATOM * sizeof(float));
      (*rotinf)[i].rot[j].Y = memY;
      memZ =  _mol_malloc (MAX_ATOM * sizeof(float));
      (*rotinf)[i].rot[j].Z = memZ;

    }
  }
}
//typedef struct rot_info * ROTPTR;

void init_rotinf(struct atomgrp *ag, int num_of_res_interface, int *res_list_interface, char *lib_contents, struct rot_info **rotinf){
  int i,j;
  int res_num;

  (*rotinf) = (struct rot_info *) _mol_malloc(num_of_res_interface*sizeof(struct rot_info));

  for(i=0;i<num_of_res_interface;i++){

    (*rotinf)[i].rot = (struct rot_list *)_mol_malloc(MAX_ROT * sizeof(struct rot_list));
   
    float* memX;
    float* memY;
    float* memZ;

    for(j = 0 ; j < MAX_ROT ; j++){

      memX =  _mol_malloc (MAX_ATOM * sizeof(float));
      (*rotinf)[i].rot[j].X = memX;
      memY =  _mol_malloc (MAX_ATOM * sizeof(float));
      (*rotinf)[i].rot[j].Y = memY;
      memZ =  _mol_malloc (MAX_ATOM * sizeof(float));
      (*rotinf)[i].rot[j].Z = memZ;

    }

    res_num = res_list_interface[i];
    (*rotinf)[i].res_num = res_num;

   get_rotamer_coordinates(ag, res_num, lib_contents, &((*rotinf)[i]), 0);

  }
}
