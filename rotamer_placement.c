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
#include _MOL_INCLUDE_
#include "main.h"
#include "parameters.h"
#include "rotamer_placement.h"
#include "graph.h"

void ifres_list_malloc( struct ifres_list** res_list )
{
        *res_list = (struct ifres_list*) _mol_malloc (sizeof (struct ifres_list));
        (*res_list)->ifres_num = (int*)  _mol_malloc ( memK * MAX_RES * sizeof(int) ) ;         // NEW ---------
}

//

void ifres_Clus_malloc( struct ifres_list** res_list , int nres_IF )
{
        *res_list = (struct ifres_list*) _mol_malloc (sizeof (struct ifres_list));
        (*res_list)->ifres_num = (int*) _mol_malloc ( nres_IF * sizeof(int) );          // NEW ---------
}

//

void Free_ifres_list( struct ifres_list** res_list )
{
	free((*res_list)->ifres_num);
	free(*res_list);
}


//

void array_sort( int* array , int size )
{
        int g,r,c;
        for ( r = 0; r < size - 1; r++){
             for ( g = r+1; g < size ;g++){
                 if ( array[r] > array[g] ){
                        c = array[r];
                        array[r] = array[g];
                        array[g] = c;
                 }
             }
        }
}

//

static void IFres_refinement( struct atomgrp* ag , struct ifres_list* res_list_in , struct ifres_list* res_list_out)
{
	int remove ; // if 1 remove the residue / OW keep it
	int resi , resj;
        int* list1 = malloc( 20 * sizeof(int));
        int* list2 = malloc( 20 * sizeof(int));
        int list1_size;
        int list2_size;	
	double rc = 5 ;
        double inter_eng;

	res_list_out->num_of_ifres = res_list_in->num_of_ifres ;
	int ind = 0 ;
	int restype;

	// RR RL LR LL :: all together

	for(int i = 0 ; i < res_list_in->num_of_ifres ; i++){
           resi = res_list_in->ifres_num[i];
	   restype = ag->res_type[resi];
	   if ( (restype != 1) && (restype != 7) && (restype != 11)  ){
		remove = -1;
		get_list( ag, resi , list1 , &list1_size);
		for(int j = 0 ; j < res_list_in->num_of_ifres ; j++){
			resj = res_list_in->ifres_num[j];
			get_list( ag, resj , list2 , &list2_size);

                        if( (i != j) && (abs(resi-resj) != 1) ){
				vdweng_inter( ag , rc, list1, list1_size, list2, list2_size, &inter_eng);
				if( abs(inter_eng) > vdw_threshold )	{  remove = 0; break; }
				else       remove = 1;
			}
		}
		if(remove == 0) {
			res_list_out->ifres_num[ind] = res_list_in->ifres_num[i];
			ind++;
		}
	   }
	}

	res_list_out->num_of_ifres = ind;

	free(list1);
	free(list2);
}


static void name_to_ID(char* res_name, double phi, double psi, int* resi_ID, int* phi_psi_ID)
{
        int res_ID;
        double phi1;
        double psi1;
        double ID_pp;

// Residue Name Converting to Residue ID

        if(strcmp(res_name,"ARG") == 0)
                res_ID=0;
        else if (strcmp(res_name,"ASN") == 0)
                res_ID=1;
        else if (strcmp(res_name,"ASP") == 0)
                res_ID=2;
        else if (strcmp(res_name,"CYS") == 0)
                res_ID=3;
        else if (strcmp(res_name,"GLN") == 0)
                res_ID=4;
        else if (strcmp(res_name,"GLU") == 0)
                res_ID=5;
        else if (strcmp(res_name,"HIS") == 0)
                res_ID=6;
	else if (strcmp(res_name,"HSD") == 0)
                res_ID=6;
        else if (strcmp(res_name,"ILE") == 0)
                res_ID=7;
        else if (strcmp(res_name,"LEU") == 0)
                res_ID=8;
        else if (strcmp(res_name,"LYS") == 0)
                res_ID=9;
        else if (strcmp(res_name,"MET") == 0)
                res_ID=10;
        else if (strcmp(res_name,"PHE") == 0)
                res_ID=11;
        else if (strcmp(res_name,"PRO") == 0)
                res_ID=12;
        else if (strcmp(res_name,"SER") == 0)
                res_ID=13;
        else if (strcmp(res_name,"THR") == 0)
                res_ID=14;
        else if (strcmp(res_name,"TRP") == 0)
                res_ID=15;
        else if (strcmp(res_name,"TYR") == 0)
                res_ID=16;
        else if (strcmp(res_name,"VAL") == 0)
                res_ID=17;
        else if (strcmp(res_name,"ALA") == 0)
                res_ID=18;
        else if (strcmp(res_name,"GLY") == 0)
                res_ID=19;
        else{
	  //printf("%s\n", "ERROR: Residue Name Is Not In Library Format.");
  	        printf("ERROR: Residue Name: %s,  Is Not In Library Format.",res_name);
                exit(EXIT_FAILURE);
        }

        resi_ID[0]=res_ID;


// Phi and Psi converting to Phi_Psi_ID

        // Rounding Phi & Psi into Library format

        if ( (phi > 180) | (phi < -180) | (psi > 180) | (psi < -180) ) {
                printf("%s\n", "ERROR: Phi or Psi angle should be in range of [-180,180].");
                exit(EXIT_FAILURE);
        }

        phi = phi/10;
        psi = psi/10;
        phi = roundf(phi);
        psi = roundf(psi);
        phi = phi * 10;
        psi = psi * 10;

        // Finding the Index

        phi1=phi/10;
        psi1=psi/10;

        phi1=phi1+18;
        psi1=psi1+18;

        ID_pp = psi1 + (37*phi1);
        phi_psi_ID[0] = (int)(ID_pp);
}

void IFres_adjustment(struct atomgrp* RecLigAg, struct ifres_list* res_list, int* cluster_i, struct rot_info* Major_rotinf, int nIFres, struct rot_info* Minor_rotinf)
{
        struct ifres_list* res_list_in;
        ifres_list_malloc( &res_list_in ) ;
        res_list_in->num_of_ifres = nIFres ;
        for(int i = 0 ; i < nIFres ; i++){
                res_list_in->ifres_num[i] = Major_rotinf[i].res_num; //// Doubt
                Minor_rotinf[i].res_num = Major_rotinf[i].res_num;
        }

        IFres_refinement( RecLigAg , res_list_in , res_list );
	// Considering full-protein model
	 
        int nresIF = res_list->num_of_ifres;

        for(int i = 0 ; i < nresIF ; i++){
                for(int j = 0 ; j < nIFres ; j++){
                        if(res_list->ifres_num[i] == res_list_in->ifres_num[j]){
                                Minor_rotinf[i].res_num = Major_rotinf[j].res_num;
                                Minor_rotinf[i].natoms = Major_rotinf[j].natoms;				
                                Minor_rotinf[i].nrotamers = Major_rotinf[j].nrotamers;
                                for(int k = 0 ; k <= Minor_rotinf[i].nrotamers ; k++){
                                        Minor_rotinf[i].rot[k].probability = Major_rotinf[j].rot[k].probability;
                                        for(int l = 0 ; l < Minor_rotinf[i].natoms ; l++){
                                                Minor_rotinf[i].rot[k].X[l] = Major_rotinf[j].rot[k].X[l];
                                                Minor_rotinf[i].rot[k].Y[l] = Major_rotinf[j].rot[k].Y[l];
                                                Minor_rotinf[i].rot[k].Z[l] = Major_rotinf[j].rot[k].Z[l];
                                        }
                                }
                                break;
                        }
                }
        }

   if(nresIF != 0) { 
      
        Free_ifres_list( &res_list_in );
        char resi_na[4];
        char* resi_name ;
        int resi_ID , ID2;
        int heavy_res_count = 0;
        int med_res_count = 0;
	
	for(int i = 0 ; i < res_list->num_of_ifres ; i++){
                for(int j = 0; j < 3 ; j++)
                        resi_na[j]=RecLigAg->idres[res_list->ifres_num[i]][j+10];
                resi_na[3]='\0';
                resi_name = resi_na ;
                name_to_ID(resi_name,0,0,&resi_ID,&ID2);
                if( (resi_ID == 0)||(resi_ID == 9) )
                        heavy_res_count++;
                else if ( (resi_ID == 4)||(resi_ID == 5)||(resi_ID == 10)||(resi_ID == 1) )
                        med_res_count++;
        }
	
        // For big graphs do clustering

        if( ( nresIF > 8 ) || ( heavy_res_count > 1) || ( med_res_count > 2 ) || ((heavy_res_count + med_res_count) > 1 ) )
                cluster_i[0] = 1;
        else
                cluster_i[0] = 0;

   }
   else
	printf("No Interface Residue Detected\n");
}

//________________________________________________________________________

// rotamer_placement.c

//________________________________________________________________________

void rot_info_malloc(struct rot_info** rotinf){

        *rotinf = (struct rot_info*) _mol_malloc (sizeof (struct rot_info));
        (*rotinf)->rot = (struct rot_list*) _mol_malloc ( MAX_ROT  * sizeof (struct rot_list));

        for(int i = 0 ; i < MAX_ROT ; i++){

                (*rotinf)->rot[i].X = (float*) _mol_malloc( MAX_ATOM * sizeof(float));
                (*rotinf)->rot[i].Y = (float*) _mol_malloc( MAX_ATOM * sizeof(float));
                (*rotinf)->rot[i].Z = (float*) _mol_malloc( MAX_ATOM * sizeof(float));
        }
}

//

void Free_rotinf(struct rot_info** rotinf)
{
        for(int i = 0 ; i < MAX_ROT ; i++){
                free((*rotinf)->rot[i].X);
                free((*rotinf)->rot[i].Y);
                free((*rotinf)->rot[i].Z);
        }
        free((*rotinf)->rot);
        free(*rotinf);

}

//===================================================================================================

int load_file_to_memory(const char *lib_filename, char **lib_content)
{
        unsigned int size = 0;
        FILE *f = fopen(lib_filename, "rb");
        if (f == NULL)
        {
                *lib_content = NULL;
                return -1; // -1 means file opening fail
        }
        fseek(f, 0, SEEK_END);
        size = ftell(f);
        fseek(f, 0, SEEK_SET);
        *lib_content = (char *)malloc(size+1);

        if (size != fread(*lib_content, sizeof(char), size, f))
        {
                free(*lib_content);
                return -2; // -2 means file reading fail
        }
        fclose(f);
        (*lib_content)[size] = 0;
        return size;
}

//===================================================================================================================================
// BLOCK II 
//===================================================================================================================================

static void get_dunbrack_rotamer(char* lib_content, int resi_ID, int phi_psi_ID, struct dunbrack_rotamer* rotamer)
{

	// Initialization

        int number_rot[20]={81,18,9,3,36,27,9,9,9,81,27,6,2,3,3,9,6,3,0,0};

	#define num_pp 1369 // for array : phi_psi

        int first_char_ind;
	int line_number;
	int rot_counter;

	char chi1ID[4];
        char chi2ID[4];
        char chi3ID[4];
        char chi4ID[4];
	char prob[8];

        int16_t chi1ID_int;
        int16_t chi2ID_int;
        int16_t chi3ID_int;
        int16_t chi4ID_int;

	double prob_fl;

	rotamer->rotamer_number = nrotCoef * number_rot[resi_ID] ;

	// Extraxting Chi's and Probabilities from Library	
	// Filling the Rotamer Struct

        for(rot_counter=0 ; rot_counter <= rotamer->rotamer_number ; rot_counter++){

	         if(resi_ID==0)
        	         line_number=phi_psi_ID * rotamer->rotamer_number + (rot_counter+1);
                 else{

			int sum = 0;
			for(int q=0;q<resi_ID;q++)
				sum = sum + num_pp * nrotCoef * number_rot[q];
			
			line_number= sum + phi_psi_ID * rotamer->rotamer_number + (rot_counter+1);	
		}

                first_char_ind = 37*(line_number-1);

		 for(int i=0;i<4;i++){

                        chi1ID[i]= lib_content[first_char_ind+8+i];
                        chi2ID[i]= lib_content[first_char_ind+13+i];
                        chi3ID[i]= lib_content[first_char_ind+18+i];
                        chi4ID[i]= lib_content[first_char_ind+23+i];
		}

                 for(int i=0;i<8;i++)
                        prob[i]=lib_content[first_char_ind+28+i];
		
                 chi1ID_int=atoi(chi1ID);
                 chi2ID_int=atoi(chi2ID);
                 chi3ID_int=atoi(chi3ID);
                 chi4ID_int=atoi(chi4ID);
                 prob_fl   =atof(prob);		

                 rotamer->chi1[rot_counter]=chi1ID_int;
                 rotamer->chi2[rot_counter]=chi2ID_int;
                 rotamer->chi3[rot_counter]=chi3ID_int;
                 rotamer->chi4[rot_counter]=chi4ID_int;
                 rotamer->probability[rot_counter]=prob_fl;		
       }
}


//===================================================================================================
// BLOCK IV 
//===================================================================================================

/* vp = v1 x v2 -vector product */

static void vprod3(double* v1, double* v2, double* vp)
{
        vp[0]=v1[1]*v2[2]-v1[2]*v2[1];
        vp[1]=v1[2]*v2[0]-v1[0]*v2[2];
        vp[2]=v1[0]*v2[1]-v1[1]*v2[0];
}



//


/* finds unit vector cd ,
   known: ab - unit vector
          bc - unit vector
          sabc - sin(abc)
          f - angle bcd
          t - torsion abcd
   output:
          sbcd - sin(f)
          cd - unit vector */

static void vcd(double* ab, double* bc, double sabc, 
             double f, double t, double* sbcd, double* cd)
{
        int i;
	double st, ct, cf, sf, e2[3], e1[3];
        vprod3(ab, bc, e2);
        for(i=0; i<3; i++)e2[i]/=sabc;
        vprod3(e2, bc, e1);
        sf=sinf(f);
        *sbcd=sf;
        st=sinf(t)*sf;
        ct=cosf(t)*sf;
        cf=cosf(f);
        for(i=0; i<3; i++)cd[i]=e1[i]*ct+e2[i]*st-bc[i]*cf;
       
}

//===================================================================================================
//===================================================================================================


/* dist = |v1-v2| -distance */


static void dist(double* v1, double* v2, double* distance)
{
	int i;
	double d, d2=0.0;
	for(i=0; i<3; i++)
	{
		d=v1[i]-v2[i];
		d2+=d*d;
	}
	*distance=sqrtf(d2);
}

//===================================================================================================
//===================================================================================================

/* an=v1 v2 v3 angle */

static void ang(double* v1, double* v2, double* v3, double* an)
{
	int i;
        double scal=0.0, d21=0.0, d23=0.0, s21, s23, ca;
	for(i=0; i<3; i++)
	{
		s21=v1[i]-v2[i];
		s23=v3[i]-v2[i];
		scal+=s21*s23;
                d21+=s21*s21;
                d23+=s23*s23;
	}
	ca=scal/sqrtf(d21)/sqrtf(d23);
	*an=acosf(ca);
}

//===================================================================================================
//===================================================================================================

/* toran= v1 v2 v3 v4 torsion angle */

static void torang(double* v1, double* v2, double* v3, double* v4, double* toran)
{
	double d23, dr, dy=0.0, dx=0.0;
        double e12[3], e34[3], e23[3], r[3], q[3];
        int i;
        for(i=0; i<3; i++)
	{
		e12[i]=v2[i]-v1[i];
                e34[i]=v4[i]-v3[i];
		e23[i]=v3[i]-v2[i];
	}
        d23=sqrtf(e23[0]*e23[0]+e23[1]*e23[1]+e23[2]*e23[2]);
        for(i=0; i<3; i++)e23[i]/=d23;
	vprod3(e12, e23, r);      
        dr=sqrtf(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
	for(i=0; i<3; i++)r[i]/=dr;
        vprod3(r, e23, q);
        for(i=0; i<3; i++)
	{
		dx+=e34[i]*q[i];
		dy+=e34[i]*r[i];
	}
	*toran=atan2f(dy,dx);
}

//========================================================================================

static void fourth_point_calc(double A[3],double B[3], double C[3], double angle_3, double tor_angle, double bond_len, double* D)
{
//      double* D; // The Fourth Point

        // Calou_ations
        double an;
        double ab[3], bc[3], sabc, sbcd, f,  t, cd[3]; 

        /* initialize sabc */
        ang(&A[0],&B[0],&C[0], &an);
        sabc=sinf(an);

        /* initialize ab */
        dist(&A[0],&B[0],&an);
        for(int j=0;j<3;j++)ab[j]=(B[j]-A[j])/an;

        /* initialize bc */
        dist(&B[0],&C[0],&an);
        for(int j=0;j<3;j++)bc[j]=(C[j]-B[j])/an;

        //  f t and length an -> f: angle between three atoms * t: torsional angle * an: bond length
        f=angle_3;
        t=tor_angle;

        double coef = 3.14159/180;
        f=f*coef;
        t=t*coef;
        an=bond_len;

        /* find unit vector cd */
        vcd(ab,bc,sabc,f,t,&sbcd,cd);

        /* find vector d */
        for(int j=0;j<3;j++)
                D[j]=C[j]+an*cd[j];

}

//===============================================================================
// BLOCK V 
//===============================================================================

static void rotamers_database(int resi_ID, double three_atoms[][7])
{

        for (int j = 0 ; j < 18 ; j++ )
               for (int k = 0 ; k < 7 ; k++)
                       three_atoms[j][k] = 0 ;


        // Filling three_atoms
        switch ( resi_ID ) {
        case 0: // ARG
              {
                double vec1[7]={1,3,4, 1.52,114.511,1,0       }; //      CG
                double vec2[7]={3,4,5, 1.52,109.797,2,0       }; //      CD
                double vec3[7]={4,5,6, 1.46,111.283,3,0       }; //      NE
                double vec4[7]={5,6,7, 1   ,115.250,4,-178.656}; //      HE
                double vec5[7]={5,6,7, 1.34,122.199,4,0};        //      CZ
                double vec6[7]={6,7,9, 1.33,120.238,0,-0.020  }; //      NH1
                double vec7[7]={7,9,10,1,   119.766,0,-0.254  }; //      HH11
                double vec8[7]={7,9,10,1,   119.766,0,180.0  };  //      HH12
                double vec9[7]={6,7,9, 1.33,119.261,0,-179.866}; //      NH2
                double vec10[7]={7,9,13,1,   119.766,0,-0.254  };//      HH21
                double vec11[7]={7,9,13,1,   119.766,0,180.0   };//      HH22

            for(int i = 0 ; i < 7 ; i++){
                three_atoms[5][i]=vec1[i];
                three_atoms[6][i]=vec2[i];
                three_atoms[7][i]=vec3[i];
                three_atoms[8][i]=vec4[i];
                three_atoms[9][i]=vec5[i];
                three_atoms[10][i]=vec6[i];
                three_atoms[11][i]=vec7[i];
                three_atoms[12][i]=vec8[i];
                three_atoms[13][i]=vec9[i];
                three_atoms[14][i]=vec10[i];
                three_atoms[15][i]=vec11[i];
            }
                }
                break;
        case 1: // ASN
                {
                double vec12[7]={1,3,4,  1.52,112.770,1,0       };
                double vec13[7]={3,4,5,  1.22,121.567,2,0       };
                double vec14[7]={3,4,5,  1.34,115.272,2,-181.419};
                double vec15[7]={4,5,7,  1   ,121.572,0,-2.978  };
                double vec16[7]={4,5,7,  1   ,117.215,0,177.740 };

                for(int i = 0 ; i < 7 ; i++){
                three_atoms[5][i]=vec12[i];
                three_atoms[6][i]=vec13[i];
                three_atoms[7][i]=vec14[i];
                three_atoms[8][i]=vec15[i];
                three_atoms[9][i]=vec16[i];
                }
                }
                break;
        case 2: // ASP
                {
                double vec17[7]={1,3,4,  1.52,112.770,1,0       };
                double vec18[7]={3,4,5,  1.24,118.005,2,0       };
                double vec19[7]={3,4,5,  1.24,118.814,2,-180.230};

            for(int i = 0 ; i < 7 ; i++){
                three_atoms[5][i]=vec17[i];
                three_atoms[6][i]=vec18[i];
                three_atoms[7][i]=vec19[i];
                }
                }
                break;
        case 3: // CYS
                {
                double vec20[7]={1,3,4,  1.77,112.771,1,0       };
            for(int i = 0 ; i < 7 ; i++){
                three_atoms[5][i]=vec20[i];
                }
                }
                break;
        case 4: // GLN
                {
                double vec21[7]={1,3,4,  1.51,112.463,1,0       };
                double vec22[7]={3,4,5,  1.51,110.311,2,0       };
                double vec23[7]={4,5,6,  1.22,120.483,3,0       };
                double vec24[7]={4,5,6,  1.34,116.468,3,-179.461};
                double vec25[7]={5,6,8,  1   ,117.215,0,176.588 };
                double vec26[7]={5,6,8, 1   ,127.272,0,-4.40   };

            for(int i = 0 ; i < 7 ; i++){
                three_atoms[5][i]= vec21[i];
                three_atoms[6][i]=vec22[i];
                three_atoms[7][i]=vec23[i];
                three_atoms[8][i]=vec24[i];
                three_atoms[9][i]=vec25[i];
                three_atoms[10][i]=vec26[i];
                }
                }
                break;
        case 5: // GLU
                {
                double vec27[7]={1,3,4, 1.51, 112.215,1,0       };
                double vec28[7]={3,4,5, 1.50, 108.585,2,0       };
                double vec29[7]={4,5,6, 1.24, 117.213,3,0       };
                double vec30[7]={4,5,6, 1.24, 117.970,3,-179.327};

            for(int i = 0 ; i < 7 ; i++){
                three_atoms[5][i]=vec27[i];
                three_atoms[6][i]=vec28[i];
                three_atoms[7][i]=vec29[i];
                three_atoms[8][i]=vec30[i];
                }
                }
                break;
        case 6: // HIS
                {
                double vec31[7]={1,3,4, 1.48, 110.670,1,0       };
                double vec32[7]={3,4,5, 1.34, 121.712,2,0       };
                double vec33[7]={4,5,6, 0.98, 126.584,0,0.851   };
                double vec34[7]={3,4,5, 1.38, 128.191,2,-182.868};
                double vec35[7]={4,5,8, 1.36, 105.126,0,176.831 };
                double vec36[7]={5,8,9,1.36, 109.299,0,0.599   };

            for(int i = 0 ; i < 7 ; i++){
                three_atoms[5][i]=vec31[i];
                three_atoms[6][i]=vec32[i];
                three_atoms[7][i]=vec33[i];
                three_atoms[8][i]=vec34[i];
                three_atoms[9][i]=vec35[i];
                three_atoms[10][i]=vec36[i];
                }
            }
                break;
        case 7: // ILE
                {

                double vec37[7]={1,3,4, 1.52, 112.580,1,-237.102};
                double vec38[7]={1,3,4, 1.52, 108.824,1,0       };
                double vec39[7]={3,4,6, 1.52, 110.332,2,0       };

            for(int i = 0 ; i < 7 ; i++){
                three_atoms[5][i]=vec37[i];
                three_atoms[6][i]=vec38[i];
                three_atoms[7][i]=vec39[i];
                }
            }
                break;
        case 8: // LEU
                {
                double vec40[7]={1,3,4, 1.50, 112.803,1,0       };
                double vec41[7]={3,4,5, 1.52, 108.329,2,0       };
                double vec42[7]={3,4,5, 1.52, 108.249,2,118.490 };

            for(int i = 0 ; i < 7 ; i++){
                three_atoms[5][i]=vec40[i];
                three_atoms[6][i]=vec41[i];
                three_atoms[7][i]=vec42[i];
                }
            }
                break;
        case 9: // LYS
                {
                double vec43[7]={1,3,4, 1.52, 109.109,1,0       };
                double vec44[7]={3,4,5, 1.51, 112.251,2,0       };
                double vec45[7]={4,5,6, 1.53, 104.178,3,0       };
                double vec46[7]={5,6,7, 1.47, 108.427,4,0       };
                double vec47[7]={6,7,8, 1.04, 109.440,0,59.819  };
                double vec48[7]={6,7,8,1.04, 109.486,0,-60.166 };
                double vec49[7]={6,7,8,1.04, 109.317,0,179.757 };

            for(int i = 0 ; i < 7 ; i++){
                three_atoms[5][i]=vec43[i];
                three_atoms[6][i]=vec44[i];
                three_atoms[7][i]=vec45[i];
                three_atoms[8][i]=vec46[i];
                three_atoms[9][i]=vec47[i];
                three_atoms[10][i]=vec48[i];
                three_atoms[11][i]=vec49[i];
                }
            }
                break;
        case 10: // MET
                {
                double vec50[7]={1,3,4, 1.53, 113.788,1,0       };
                double vec51[7]={3,4,5, 1.80, 113.005,2,0       };
                double vec52[7]={4,5,6, 1.80, 100.296,3,0       };

            for(int i = 0 ; i < 7 ; i++){
                three_atoms[5][i]=vec50[i];
                three_atoms[6][i]=vec51[i];
                three_atoms[7][i]=vec52[i];
                }
            }
                break;
        case 11: // PHE
                {
                double vec53[7]={1,3,4, 1.48, 112.545,1,0       };
                double vec54[7]={3,4,5, 1.40, 119.767,2,0       };
                double vec55[7]={3,4,5, 1.40, 118.886,2,180.215 };
                double vec56[7]={4,5,6, 1.40, 119.067,0,-179.839};
                double vec57[7]={4,5,7, 1.40, 119.642,0,179.860 };
                double vec58[7]={5,6,8,1.40, 119.743,0,-0.185  };

            for(int i = 0 ; i < 7 ; i++){
                three_atoms[5][i]=vec53[i];
                three_atoms[6][i]=vec54[i];
                three_atoms[7][i]=vec55[i];
                three_atoms[8][i]=vec56[i];
                three_atoms[9][i]=vec57[i];
                three_atoms[10][i]=vec58[i];
                }
            }
                break;
        case 12: // PRO
                {
                double vec59[7]={3,4,5, 1.52, 107.286,2,0       };
                double vec60[7]={1,3,4, 1.52, 103.830,1,0       };

            for(int i = 0 ; i < 7 ; i++){
                three_atoms[2][i]=vec59[i];
                three_atoms[5][i]=vec60[i];
                }
            }
                break;
        case 13: // SER
                {
                double vec61[7]={1,3,4, 1.40, 104.688,1,0       };
                double vec62[7]={3,4,5, 0.95, 106.284,0,173.552 };

            for(int i = 0 ; i < 7 ; i++){
                three_atoms[5][i]=vec61[i];
                three_atoms[6][i]=vec62[i];
                }
            }
                break;
        case 14: // THR
                {
                double vec63[7]={1,3,4, 1.40, 106.364,1,0       };
                double vec64[7]={3,4,5, 0.95, 111.908,0,175.415 };
                double vec65[7]={1,3,4, 1.50, 106.476,1,-244.004};

            for(int i = 0 ; i < 7 ; i++){
                three_atoms[5][i]=vec63[i];
                three_atoms[6][i]=vec64[i];
                three_atoms[7][i]=vec65[i];
                }
            }
                break;
        case 15: // TRP
                {
                double vec66[7]={1,3,4, 1.49, 111.433,1,0       };
                double vec67[7]={3,4,5, 1.43, 126.924,2,-178.746};
                double vec68[7]={4,5,6, 1.43, 106.953,0,178.704 };
                double vec69[7]={4,5,6, 1.40, 133.050,0,-1.270  };
                double vec70[7]={3,4,5, 1.39, 127.154,2,0       };
                double vec71[7]={4,5,9,1.36, 110.569,0,-178.734};
                double vec72[7]={5,9,10,0.98,125.463,0,179.690 };
                double vec73[7]={5,6,7,1.40, 121.080,0,-179.773};
                double vec74[7]={5,6,8,1.40, 118.695,0,179.846 };
                double vec75[7]={6,8,13,1.41,120.903,0,-0.019  };

            for(int i = 0 ; i < 7 ; i++){
                three_atoms[5][i]=vec66[i];
                three_atoms[6][i]=vec67[i];
                three_atoms[7][i]=vec68[i];
                three_atoms[8][i]=vec69[i];
                three_atoms[9][i]=vec70[i];
                three_atoms[10][i]=vec71[i];
                three_atoms[11][i]=vec72[i];
                three_atoms[12][i]=vec73[i];
                three_atoms[13][i]=vec74[i];
                three_atoms[14][i]=vec75[i];
                }
            }
                break;
        case 16: // TYR
                {
                double vec76[7]={1,3,4, 1.51, 111.945,1,0       };
                double vec77[7]={3,4,5, 1.40, 120.042,2,0       };
                double vec78[7]={4,5,6, 1.40, 119.734,0,-179.763};
                double vec79[7]={3,4,5, 1.40, 119.126,2,-180.494};
                double vec80[7]={4,5,8, 1.40, 119.722,0,179.860 };
                double vec81[7]={5,6,7,1.40, 118.858,0,0.150   };
                double vec82[7]={6,7,10,1.37,120.093,0,-179.77 };
                double vec83[7]={7,10,11,0.95,109.638,0,176.521};

            for(int i = 0 ; i < 7 ; i++){
                three_atoms[5][i]=vec76[i];
                three_atoms[6][i]=vec77[i];
                three_atoms[7][i]=vec78[i];
                three_atoms[8][i]=vec79[i];
                three_atoms[9][i]=vec80[i];
                three_atoms[10][i]=vec81[i];
                three_atoms[11][i]=vec82[i];
                three_atoms[12][i]=vec83[i];
                }
            }
                break;
        case 17: // VAL
                {
                double vec84[7]={1,3,4, 1.51, 110.842 ,1,0      };
                double vec85[7]={1,3,4, 1.53, 109.728,1,-238.667};

            for(int i = 0 ; i < 7 ; i++){
                three_atoms[5][i]=vec84[i];
                three_atoms[6][i]=vec85[i];
                }
            }
                break;

        default:
                printf("ERROR on filling the three_atoms information\n");
                break;
        }

}


//===================================================================================================
// BLOCK VII  
//===================================================================================================

static void phi_psi_extraction(struct atomgrp* ag , double* phi, double* psi)
{
        double pi = 3.14159265;

        int resi_num = ag->nres ;

        double C0[3];
        double N[3];
        double CA[3];
        double C[3];
        double N1[3];

        double myphi;
        double mypsi;

        int num_prev_atom;
        int num_last_atom;

        for(int i = 0 ; i < resi_num ; i++){

                if( i == 0 )
                        num_prev_atom = 0;
                else
                        num_prev_atom = ag->iares[i];                       // Last Previous Residue Atom Number

                num_last_atom = ag->iares[i+1];                             // Last Current Residue Atom Number


                if( i != 0 ){

                        C0[0]= ag->atoms[num_prev_atom - 2 ].X;
                        C0[1]= ag->atoms[num_prev_atom - 2 ].Y;
                        C0[2]= ag->atoms[num_prev_atom - 2 ].Z;
                }
                        N[0] = ag->atoms[num_prev_atom ].X;
                        N[1] = ag->atoms[num_prev_atom ].Y;
                        N[2] = ag->atoms[num_prev_atom ].Z;

                        CA[0]= ag->atoms[num_prev_atom + 2].X;
                        CA[1]= ag->atoms[num_prev_atom + 2].Y;
                        CA[2]= ag->atoms[num_prev_atom + 2].Z;

                        C[0]= ag->atoms[num_last_atom - 2].X;
                        C[1]= ag->atoms[num_last_atom - 2].Y;
                        C[2]= ag->atoms[num_last_atom - 2].Z;

                if(i != resi_num -1){

                        N1[0] = ag->atoms[num_last_atom ].X;
                        N1[1] = ag->atoms[num_last_atom ].Y;
                        N1[2] = ag->atoms[num_last_atom ].Z;
                }

                if(i == 0){

                        phi[i] = 0;
                        torang( N , CA , C , N1 , &mypsi);
                        mypsi = mypsi * 180 / pi ;
                        psi[i] = mypsi ;
                }
                else if (i == resi_num -1){

                        torang( C0 , N , CA , C , &myphi);
                        myphi = myphi * 180 / pi ;
                        phi[i] = myphi ;
                        psi[i]= 0;
                }
                else {
                        torang( C0 , N , CA , C , &myphi);
                        myphi = myphi * 180 / pi ;
                        phi[i] = myphi ;

                        torang( N , CA , C , N1 , &mypsi);
                        mypsi = mypsi * 180 / pi ;
                        psi[i] = mypsi ;
                }
        }
}

//=====================================

void get_rotamer_coordinates( struct atomgrp* ag,  int resi_num, char* lib_contents, struct rot_info* rotinf, int number_rotamers )
{

// Phi-Psi Extraction

        double* myphi;
        double* mypsi;
        myphi = malloc(MAX_RES*sizeof(double));
        mypsi = malloc(MAX_RES*sizeof(double));

        phi_psi_extraction(ag , myphi, mypsi);

// Converting Phi_Psi into Library format

        double phi = myphi[resi_num];
        double psi = mypsi[resi_num];

        phi = phi / 10;
        psi = psi / 10;
        phi = roundf(phi);
        psi = roundf(psi);
        phi = phi * 10;
        psi = psi * 10;

// Finding the Index: Coding the phi_psi

        double phi1=phi/10;
        double psi1=psi/10;

        phi1=phi1+18;
        psi1=psi1+18;

        double ID_pp = psi1 + (37*phi1);
        int phi_psi_ID = (int)(ID_pp);

        free(myphi);
        free(mypsi);

// Finding Residue Name

        char resi_na1[4];

       for(int j = 0; j < 3 ; j++)
                resi_na1[j]=ag->idres[resi_num][j+10];
       resi_na1[3]='\0';

       char* resid_name = resi_na1 ;

// Reading Data from Atomgroup Struct

       int* Res_List;
       Res_List = malloc( ag->nres * sizeof(int));

       char resi_na[4];
       int resi_id;
       int ID3;
       int num_of_resi = ag->nres ;

       for(int i = 0; i < num_of_resi ; i++){
                for(int k = 0; k < 3 ; k++)
                        resi_na[k]=ag->idres[i][k+10];
                resi_na[3]='\0';

                name_to_ID(resi_na, 0 , 0 , &resi_id, &ID3);

                Res_List[i] = resi_id;        // Res_List[0]=0; The first Resiue is ARG
        }

// Start Placing the Atoms

        int i,j,k; 

        int num_atoms_res[20]={17,11,9,7,12,10,12,9,9,13,9,12,7,8,9,16,14,8,6,5};  // num_atoms_res(0)=17 -> ARG
        int number_rot[20]={81,18,9,3,36,27,9,9,9,81,27,6,2,3,3,9,6,3,0,0};        // num_rotamers_res(0)=81 -> ARG


        int resi_ID;            // LYS: 9
        int resi_natoms;        // 13 number of LYS atoms
        int resi_nrots;         // 81 number of rotamers of ARG
        int ID2;

        name_to_ID(resid_name,0,0,&resi_ID,&ID2); // resi_ID : convert LYS into 9
        resi_natoms=num_atoms_res[resi_ID]; // Number of atoms of LYS: 13
        resi_nrots= nrotCoef * number_rot[resi_ID];

// Three_Atoms: for finding atom orders and residue angles and bond_length information

        // Fixed atoms: 1,2,3,4, and last two

	double three_atoms[18][7];
             

// Extracting Library Information: Chi's and Probability

        struct dunbrack_rotamer rotamer;

        rotamer.chi1 = (int16_t*) _mol_malloc( MAX_ROT * sizeof(int16_t));     
        rotamer.chi2 = (int16_t*) _mol_malloc( MAX_ROT * sizeof(int16_t));    
        rotamer.chi3 = (int16_t*) _mol_malloc( MAX_ROT * sizeof(int16_t));
        rotamer.chi4 = (int16_t*) _mol_malloc( MAX_ROT * sizeof(int16_t));
        rotamer.probability = (double*) _mol_malloc( MAX_ROT * sizeof(double));

        if( (resi_ID != 18) & (resi_ID != 19) & (resi_ID != 3) )
                get_dunbrack_rotamer(lib_contents, resi_ID, phi_psi_ID, &rotamer);

// Memory Allocation and copying coordinate from Atom Group

        double XX[resi_nrots][resi_natoms];       //XX[2][3] -> x of 4th atom of 3rd rotamer from given residue
        double YY[resi_nrots][resi_natoms];
        double ZZ[resi_nrots][resi_natoms];

        int num_prev_atom = ag->iares[resi_num];

        for(i=0;i< resi_nrots;i++){
                for(j=0;j<resi_natoms;j++){
                        XX[i][j]= ag->atoms[num_prev_atom + j].X;
                        YY[i][j]= ag->atoms[num_prev_atom + j].Y;
                        ZZ[i][j]= ag->atoms[num_prev_atom + j].Z;
                }
        }

// FOR Non-Bonund Structure

        for(j = 0 ; j < resi_natoms ; j++){

                rotinf->rot[0].X[j] = XX[0][j];
                rotinf->rot[0].Y[j] = YY[0][j];
                rotinf->rot[0].Z[j] = ZZ[0][j];

        }
// Fillinf rotinf

        rotinf->natoms    = num_atoms_res[resi_ID];
        rotinf->nrotamers = nrotCoef * number_rot[resi_ID];

// Main Loop to Placing Atoms in Rotamers

        int   ind1, ind2, ind3;
        double A[3];
        double B[3];
        double C[3];
        double D[3];
        double bond_len;
        double angle3;
        int chi_ind;
        double delta;
        int chi_ang_ind;
        double chi_ang_ind_f;
        double chi_ang;
        double tor_ang = 0;
        
        // Naming the output files

        int rot_num_treshold;

        if(number_rotamers < 0)
                rot_num_treshold = 0 ;
        else if(number_rotamers == 0)
                rot_num_treshold = resi_nrots ;
        else if(number_rotamers > resi_nrots)
                rot_num_treshold = resi_nrots;
        else
                rot_num_treshold = number_rotamers;

        if( (resi_ID != 18) & (resi_ID != 19) & (resi_ID != 3) )
                rotamers_database(resi_ID, three_atoms);

      for(i = 0 ; i < rot_num_treshold ; i++){

           if( (resi_ID != 18) & (resi_ID != 19) & (resi_ID != 3)  ){

                if(resi_ID == 12){      // PRO: different from others

                       k=5;
                       int ii;
                       for(ii = 0 ; ii < 2 ; ii++){  // at (k+1)th atom

                                        ind1 = (int) three_atoms[k][0];        // First Atom Index
                                        ind2 = (int) three_atoms[k][1];        // Second Atom Index
                                        ind3 = (int) three_atoms[k][2];        // Third Atom Index

                                        A[0]= XX[i][ind1-1];
                                        A[1]= YY[i][ind1-1];
                                        A[2]= ZZ[i][ind1-1];

                                        B[0]= XX[i][ind2-1];
                                        B[1]= YY[i][ind2-1];
                                        B[2]= ZZ[i][ind2-1];

                                        C[0]= XX[i][ind3-1];
                                        C[1]= YY[i][ind3-1];
                                        C[2]= ZZ[i][ind3-1];

                                        bond_len = three_atoms[k][3];
                                        angle3=three_atoms[k][4];
                                        chi_ind = (int) three_atoms[k][5];
                                        delta = three_atoms[k][6];

                                        switch ( chi_ind  ) {
                                        case 0:{
                                                        tor_ang =  delta;
                                                }
                                                break;
                                        case 1:{
                                                        chi_ang_ind = rotamer.chi1[i];
                                                        chi_ang_ind = chi_ang_ind - 1800;
                                                        chi_ang_ind_f = (double) chi_ang_ind;
                                                        chi_ang =  chi_ang_ind_f / 10;

                                                        tor_ang = chi_ang + delta;
                                                }
                                                break;
                                        case 2:{
                                                        chi_ang_ind = rotamer.chi2[i];
                                                        chi_ang_ind = chi_ang_ind - 1800;
                                                        chi_ang_ind_f = (double) chi_ang_ind;
                                                        chi_ang =  chi_ang_ind_f / 10;

                                                        tor_ang = chi_ang + delta;
                                                }
                                                break;
                                        case 3:{
                                                        chi_ang_ind = rotamer.chi3[i];
                                                        chi_ang_ind = chi_ang_ind - 1800;
                                                        chi_ang_ind_f = (double) chi_ang_ind;
                                                        chi_ang = chi_ang_ind_f / 10;

                                                        tor_ang = chi_ang + delta;
                                                }
                                                break;
                                        case 4:{
                                                        chi_ang_ind = rotamer.chi4[i];
                                                        chi_ang_ind = chi_ang_ind - 1800;
                                                        chi_ang_ind_f = (double) chi_ang_ind;
                                                        chi_ang =  chi_ang_ind_f / 10;

                                                        tor_ang = chi_ang + delta;
                                                }
                                                break;

                                        default:
                                                printf("Error in reading the torsional angle from Struct\n");
                                                break;
                                                }

                                        fourth_point_calc( A,B,C , angle3 , tor_ang, bond_len, D);

                                        // Updating the PDB
                                        XX[i][k-1]=D[0];
                                        YY[i][k-1]=D[1];
                                        ZZ[i][k-1]=D[2];

                                        k=2;

                        }
                }

                else{

                        for(k = 5 ; k <= resi_natoms-2 ; k++){  // at (k+1)th atom

                                        ind1 = (int) three_atoms[k][0];        // First Atom Index
                                        ind2 = (int) three_atoms[k][1];        // Second Atom Index
                                        ind3 = (int) three_atoms[k][2];        // Third Atom Index

                                        A[0]= XX[i][ind1-1];
                                        A[1]= YY[i][ind1-1];
                                        A[2]= ZZ[i][ind1-1];

                                        B[0]= XX[i][ind2-1];
                                        B[1]= YY[i][ind2-1];
                                        B[2]= ZZ[i][ind2-1];

                                        C[0]= XX[i][ind3-1];
                                        C[1]= YY[i][ind3-1];
                                        C[2]= ZZ[i][ind3-1];

                                        bond_len = three_atoms[k][3];
                                        angle3=three_atoms[k][4];
                                        chi_ind = (int) three_atoms[k][5];
                                        delta = three_atoms[k][6];

                                        switch ( chi_ind  ) {
                                        case 0:{
                                                        tor_ang =  delta;
                                                }
                                                break;
                                        case 1:{
                                                        chi_ang_ind = rotamer.chi1[i];
                                                        chi_ang_ind = chi_ang_ind - 1800;
                                                        chi_ang_ind_f = (double) chi_ang_ind;
                                                        chi_ang =  chi_ang_ind_f / 10;

                                                        tor_ang = chi_ang + delta;
                                                }
                                                break;
                                        case 2:{
                                                        chi_ang_ind = rotamer.chi2[i];
                                                        chi_ang_ind = chi_ang_ind - 1800;
                                                        chi_ang_ind_f = (double) chi_ang_ind;
                                                        chi_ang =  chi_ang_ind_f / 10;

                                                        tor_ang = chi_ang + delta;
                                                }
                                                break;
                                        case 3:{
                                                        chi_ang_ind = rotamer.chi3[i];
                                                        chi_ang_ind = chi_ang_ind - 1800;
                                                        chi_ang_ind_f = (double) chi_ang_ind;
                                                        chi_ang = chi_ang_ind_f / 10;

                                                        tor_ang = chi_ang + delta;
                                                }
                                                break;
                                        case 4:{
                                                        chi_ang_ind = rotamer.chi4[i];
                                                        chi_ang_ind = chi_ang_ind - 1800;
                                                        chi_ang_ind_f = (double) chi_ang_ind;
                                                        chi_ang =  chi_ang_ind_f / 10;

                                                        tor_ang = chi_ang + delta;
                                                }
                                                break;

                                        default:
                                                printf("Error in reading the torsional angle from Struct\n");
                                                break;
                                                }

                                        fourth_point_calc( A,B,C , angle3 , tor_ang, bond_len, D);

                                        // Updating the PDB

                                        XX[i][k-1]=D[0];
                                        YY[i][k-1]=D[1];
                                        ZZ[i][k-1]=D[2];

                        }
                 }
           }

                for(j = 0 ; j < resi_natoms ; j++){

                        rotinf->rot[i+1].X[j] = XX[i][j];
                        rotinf->rot[i+1].Y[j] = YY[i][j];
                        rotinf->rot[i+1].Z[j] = ZZ[i][j];
                }

		rotinf->rot[i+1].probability = rotamer.probability[i];		

	}

        if( (resi_ID != 18) & (resi_ID != 19) & (resi_ID != 3)  )
                rotinf->rot[0].probability = rotinf->rot[1].probability + rotinf->rot[1].probability*0.2;
	else
		rotinf->rot[0].probability = 1;

        
        free( rotamer.chi1 ) ;
        free( rotamer.chi2 ) ; 
        free( rotamer.chi3 ) ; 
        free( rotamer.chi4 ) ; 
        free( rotamer.probability ) ; 
        free(Res_List);         

}

//===================================================================================================

void apply_rot_to_atmgrp( struct atomgrp* ag,  struct rot_info* rotinf, int resi_number, int rotamer_index )
{
        int num_atoms_res[20]={17,11,9,7,12,10,12,9,9,13,9,12,7,8,9,16,14,8,6,5};  // num_atoms_res(0)=17 -> ARG

// Updating the AG with the rotamer

        int num_prev_atom = ag->iares[resi_number];

// Finding Residue Name

       char resi_na[4];

       for(int j = 0; j < 3 ; j++)
                resi_na[j]=ag->idres[resi_number][j+10];
       resi_na[3]='\0';

       char* resi_name = resi_na ;

// Convert The name to ID and Finding the number of atoms of residue

       int resi_ID;
       int ID2;
       name_to_ID(resi_name,0,0,&resi_ID,&ID2);
       int resi_natoms=num_atoms_res[resi_ID];

// Updating the Atom Group

        if ( (resi_ID != 18) & (resi_ID != 19) & (resi_ID != 3)  ){
                for(int i = 4 ; i < resi_natoms - 2 ; i++){
                       ag->atoms[num_prev_atom + i].X = rotinf->rot[rotamer_index].X[i];
                       ag->atoms[num_prev_atom + i].Y = rotinf->rot[rotamer_index].Y[i];
                       ag->atoms[num_prev_atom + i].Z = rotinf->rot[rotamer_index].Z[i];
                }
        }
}
