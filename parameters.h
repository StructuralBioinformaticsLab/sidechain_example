
// parameters.h

#ifndef _PARAMETERS_H
#define _PARAMETERS_H

/* Interface Parameters */
double cutoff;
#define vdw_threshold 10.0
#define lig_rec_dist  12.0	// to mark interfacial residues

/* Energy Terms Parametes */

#define vdw_noIF_coef   0.8 	
#define vdw_inter_coef  1	
#define prob1_coef      42.5
#define hbond_coef      75.3
#define nb_rot_prob     0.1

/* Refinement Parameters */
#define rmsd_trshold   -1     // RMSD for Rotamer Compression

/* Library Switch */
#define LIB_i 0         // if 1 : extended library - if 0 : original library
int nrotCoef;

/* Controlling Parametes */
#define big_eng         10000
#define Small_Max_Node    200 
#define Big_Max_Node     2000
#define Global_Max_Node   300	// def : 200
#define max_add1          0.6 
#define max_add2          0.2

/* Fixed Parameters */
int MAXNODE ;
int MAX_RES ;
int MAX_ROT;
double LARGE ;
double MINIM ;
#define P0 	     1    // For Boltzman Energy
#define memK 	   100    // For Interface List Memory Allocation
#define MAX_LEN    200
#define MAX_ATOM    20
#define MAX_CHAR   200

/* Function Indicators */
#define DEE_i  1 // if 1 : do DEE - if 0 : do not
#define COMP_i 0 // for Rotamer_Compression
#define UNB_i  1 // if 1 : do not consider unbound rotamers

/* Files */
char* prmfile ;
char* rtffile ;

/* CYS */
int* CYS_LIST;
int CYS_L_SIZE;

/* Rotamer_Filter based on Probability */

#define  Prob_Sum 0.9
#define chi_tresh 40 

/* Setting for MWIS.c */

#define MWIS_mode 0       // if 0: mwis-v1 :: fixed step size
			  // if 1: mwis-v2 :: dynamic step size

#define myTHRESHOLD1       0.4
#define myTHRESHOLD0       0.2
#define N_mw1              20000  	
#define delta_mw1          0.3
#define lambda_mw1         0.3
#define epsilon_mw1        1.0
#define eps_lowbound_mw1   0.1

#define N_mw2		   10000
#define epsilon_mw2 	   1000	
#define eps_lowbound_mw2   100


#endif
