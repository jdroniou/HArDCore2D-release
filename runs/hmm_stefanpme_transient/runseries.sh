#!/bin/bash
#
# Execute hho file on series of meshes, and calculate outputs
#

###
# Executable file for the considered scheme
executable_name="hmm-stefanpme-transient"

###
# Directories
origin=$(pwd)
if [ ! -f ../directories.sh ]; then
  echo "directories.sh does not exist. Please read the README.txt in the parent folder."
  exit
fi
. ../directories.sh

# Options:
if [[ $1 == "help" ]]; then
	echo -e "\nExecute tests using parameters in data.sh, creates and compile latex file, and calculate rates.\n"
  exit;
fi;


if [ ! -d $outdir ]; then
  mkdir $outdir
fi
\rm -r $outdir/*

###
# LOAD DATA
# (LATER: Test that each required data exists (files, parameters...))
. data.sh

echo "boundary conditions bc=$bc"
echo -e "test case: solution=$tcsol, diffusion=$tcdiff\n"

# Create sequence of time steps
nbmesh=${#mesh[@]}
dt[1]=$dt_initial
for i in `seq 2 $nbmesh`; do
  dt[$i]=$(echo "(${dt[$(($i-1))]}/4.0)/1" | scale=0 bc -l);
done;

###
# EXECUTE FOR EACH MESH
for i in `seq 1 $nbmesh`; 
do
  # Mesh, and time step divided by 2 at each new iteration
  meshfile=$meshdir"/"${mesh[$i]}".typ2"
  timestep=${dt[$i]}
  echo -e "\n*************************\nmesh $i out of $nbmesh: ${mesh[$i]}.typ2 (time step: $timestep)"
  # Execute code
  $executable -m $meshfile --FinalTps $FinalTps --dt $timestep -b $bc -c $tcsol $tcdiff --testcaseNL $tcNL --powerPME $mPME --weight $weight --source $source --solver_type $solver_type --use_threads $use_threads
  r=$?
  if [ "$r" != "0" ]; then
    exit
  fi
  # Move outputs
  mv results.txt $outdir/results-$i.txt
  if [ -f T-solution.vtu ]; then
    mv T-solution.vtu $outdir/T-solution-$i.vtu
  fi
  mv solution.vtu $outdir/mesh"$i"_solution.vtu
  mv zeta-solution.vtu $outdir/mesh"$i"_zeta-solution.vtu
  mv exact-solution.vtu $outdir/mesh"$i"_exact-solution.vtu
  mv zeta-exact-solution.vtu $outdir/mesh"$i"_zeta-exact-solution.vtu
done;

# CREATE DATA FILE FOR LATEX
cd $outdir
echo -e "meshsize timestep Lmp1error rate H1error rate EnergyError NbEdgeDOFs AveNewton" > $errorsfile
for i in `seq 1 $nbmesh`; 
do
  MeshSize=$(awk '/MeshSize:/ {print $NF}' results-$i.txt)
  TimeStep=$(awk '/TimeStep:/ {print $NF}' results-$i.txt)
  L2error=$(awk '/L2error:/ {print $NF}' results-$i.txt)
  Lmp1error=$(awk '/Lmp1error:/ {print $NF}' results-$i.txt)
  H1error=$(awk '/H1error:/ {print $NF}' results-$i.txt)
  EnergyError=$(awk '/EnergyError:/ {print $NF}' results-$i.txt)
  NbEdgeDOFs=$(awk '/NbEdgeDOFs:/ {print $NF}' results-$i.txt)
  MeshReg=$(awk '/MeshReg:/ {print $NF}' results-$i.txt)
  MeshAnisotropy=$(awk '/MeshAnisotropy:/ {print $NF}' results-$i.txt)
  AveNewton=$(awk '/AveNewton:/ {print $NF}' results-$i.txt)
  FractionNegValues=$(awk '/FractionNegValues:/ {print $NF}' results-$i.txt)
  FractionNegMass=$(awk '/FractionNegMass:/ {print $NF}' results-$i.txt)
  MinValue=$(awk '/MinValue:/ {print $NF}' results-$i.txt)
  MaxValue=$(awk '/MaxValue:/ {print $NF}' results-$i.txt)      
  if(($i > 1)); then
    imo=$(perl -E "say $i - 1")
    OldMeshSize=$(awk '/MeshSize:/ {print $NF}' results-$imo.txt)
    OldLmp1error=$(awk '/Lmp1error:/ {print $NF}' results-$imo.txt)
    Lmp1ErrorRate=$(perl -E "say sprintf(\"%.2f\", log($OldLmp1error/$Lmp1error)/log($OldMeshSize/$MeshSize))")
    OldH1error=$(awk '/H1error:/ {print $NF}' results-$imo.txt)
    H1ErrorRate=$(perl -E "say sprintf(\"%.2f\", log($OldH1error/$H1error)/log($OldMeshSize/$MeshSize))")
    echo -e "$MeshSize $TimeStep $Lmp1error $Lmp1ErrorRate $H1error $H1ErrorRate $EnergyError $NbEdgeDOFs $AveNewton" >> $errorsfile
  else
    echo -e "$MeshSize $TimeStep $Lmp1error -- $H1error -- $EnergyError $NbEdgeDOFs $AveNewton" >> $errorsfile
  fi
done;

column -t < $errorsfile | tee $errorsfile

cd ..

