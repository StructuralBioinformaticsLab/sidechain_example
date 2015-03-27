#ifndef MCM_H
#define MCM_H

#define MAX_NUM_ENG_COMPONENTS                   15
#define MC_OUT_STEPS                             10
#define KT                                       0.591//in kcal/mol
#define ADAPTIVE_TRIAL_MOVE_STEPS                10
#define whbond                                   0.441
#define wdars                                    0.05
#define wsolnp                                   0.279
#define wsolpol                                  0
#define wvdwatr                                  0.338
#define wvdwrep                                  0.338
#define wcoul                                    0.025

struct List {
  int n;
  int* K;
};

struct sort_list{
  int num;
  double eng;
};


struct rot_info {
  int res_num;//residue number of interfacial residue
  int natoms;
  int nrotamers;
  struct rot_list *rot;
};


struct rot_list{
  float *X;
  float *Y;
  float *Z;
  float probability;
  float eng;//ln of probability
};


struct side_chain_centroids{
  int num;//the number of centroids
  double *X;
  double *Y;
  double *Z;
  double *origX;
  double *origY;
  double *origZ;
  double *GX;
  double *GY;
  double *GZ;
  double *eps;
  double *rminh;
  struct List lig_list;//a list of centroids that belong to ligand
};

struct my_par_rigid{ 
  int nbupdate; 
  struct atomgrp *ag; 
  struct agsetup *ags;
  struct acesetup *acs;
  struct prm *atomprm;
  double efacv; 
  double eface; 
  struct List lig_list;
  double *origin;
  double *center;
  double *minv;//rotation and translation vector of ligand after minimization
  double rotprob;
  double hbond;
  double dars;
  double solnp;
  double solpol;
  double vdweng_rep;
  double vdweng_atr;
  double couleng;
  double beng;
  double aeng;
  double ieng;
  double teng;
  int rotmin_res_num;
  int num_of_side_chain_atoms;
  struct side_chain_centroids *centroids;
  int tbond;// tbond is used in hydrogen minimization , if tbond = 0 we trun off teng during minimization 
}; 

#endif
