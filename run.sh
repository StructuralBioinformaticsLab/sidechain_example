
# This bash script runs our MWIS-based side-chain packing algoritm on a given protein 

if [ $# -ne 0 ]
then
	echo "Usage:"
	echo "./run.sh"
	exit 1
fi

./main 1dfj/1dfj_before_repack.pdb 1dfj/1dfj_r_l_u_nmin.psf 1dfj/1dfj_r_l_u_nmin.mol2 1dfj/1dfj_r_l_u_nmin_l.pdb prms/atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec prms/parm_new.prm prms/pdbamino_new.rtf prms/rotamer_lib.txt 1dfj/1dfj_mwis_repacked.pdb

echo "Running side-chain packing on protein with PDB ID of 1dfj" 
echo "The output file is stored in the following path:"
echo 1dfj/1dfj_mwis_repacked.pdb
