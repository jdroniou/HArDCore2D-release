#!/bin/bash
#
# Execute hho file on series of meshes, and calculate outputs
#

# Options
if [[ $1 == "help" ]]; then
    echo -e "\nExecute tests using parameters in data.sh"
    exit;
fi;

# Load data
. data.sh

# File for times
timesfile="times.dat"

# Directories
executable_name="ddr-klplate"
#make $executable_name

origin=$(pwd)
if [ ! -f ../directories.sh ]; then
  echo "directories.sh does not exist. Please read the README.txt in the parent folder."
  exit
fi
. ../directories.sh

echo "Degree                 : $k"
echo "Mesh family            : $mesh_family"
echo "Stabilization parameter: $stab_par"
echo "D/nu                   : $D/$nu"
echo "Test case: solution    : $solution"
echo "Output directory       : $outdir"

outsubdir=${outdir}/${mesh_family}_k${k}
if [ ! -d $outsubdir ]; then
    mkdir -p $outsubdir
else
    \rm -r $outsubdir/*
fi

###
# EXECUTE FOR EACH MESH
nbmesh=${#mesh[@]}
for i in `seq 1 $nbmesh`; 
do
    meshtype=$(echo ${mesh[$i]} | cut -d ':' -f 1)
    meshfile=$meshdir"/"$(echo ${mesh[$i]} | cut -d ':' -f 2)
    echo -e "------------------------------------------------------------------------------"
    echo -e "Mesh $i out of $nbmesh: $meshfile"
    echo -e "Output directory: $outsubdir"
    echo -e "------------------------------------------------------------------------------"
    # Execute code
    if($executable -m $meshfile.typ2 -k $k -s $solution -D $D -n $nu -x $stab_par -e $export_matrix); then
      # Move outputs
      mv results.txt $outsubdir/results-$i.txt
      if [ -f displacement.vtu ]; then
      	mv displacement.vtu $outsubdir/mesh"$i"_displacement.vtu
    	fi
    	if [ -f exact-displacement.vtu ]; then
      	mv exact-displacement.vtu $outsubdir/mesh"$i"_exact-displacement.vtu
      fi
    	if [ $export_matrix ]; then
      	mv A.mtx $outsubdir/A"$i".mtx
      	mv b.mtx $outsubdir/b"$i".mtx
      fi
    fi
done

# CREATE DATA FILE FOR LATEX
echo -e "Deg MeshSize NbCells NbEdges NbVertices DimXDivDiv Error Rate" > $outsubdir/$errorsfile
echo -e "TwallPlatesCore TprocPlatesCore TwallModel TprocModel TwallSolve TprocSolve" > $outsubdir/$timesfile
for i in `seq 1 $nbmesh`; 
do
    Degree=$(awk '/Degree:/ {print $NF}' $outsubdir/results-$i.txt)
    MeshSize=$(awk '/MeshSize:/ {print $NF}' $outsubdir/results-$i.txt)
    NbCells=$(awk '/NbCells:/ {print $NF}' $outsubdir/results-$i.txt)
    NbEdges=$(awk '/NbEdges:/ {print $NF}' $outsubdir/results-$i.txt)
    NbVertices=$(awk '/NbVertices:/ {print $NF}' $outsubdir/results-$i.txt)
    DimXDivDiv=$(awk '/DimXDivDiv:/ {print $NF}' $outsubdir/results-$i.txt)
    DimSystem=$(awk '/DimSystem:/ {print $NF}' $outsubdir/results-$i.txt)
    Error=$(awk '/^Error:/ {print $NF}' $outsubdir/results-$i.txt)
    TwallPlatesCore=$(awk '/TwallPlatesCore:/ {print $NF}' $outsubdir/results-$i.txt)
    TprocPlatesCore=$(awk '/TprocPlatesCore:/ {print $NF}' $outsubdir/results-$i.txt)
    TwallModel=$(awk '/TwallModel:/ {print $NF}' $outsubdir/results-$i.txt)
    TprocModel=$(awk '/TprocModel:/ {print $NF}' $outsubdir/results-$i.txt)
    TwallSolve=$(awk '/TwallSolve:/ {print $NF}' $outsubdir/results-$i.txt)
    TprocSolve=$(awk '/TprocSolve:/ {print $NF}' $outsubdir/results-$i.txt)
    echo -e "$TwallPlatesCore $TprocPlatesCore $TwallModel $TprocModel $TwallSolve $TprocSolve" >> $outsubdir/$timesfile    
    if(($i \> 1)); then
      imo=$(perl -E "say $i - 1")
      OldMeshSize=$(awk '/MeshSize:/ {print $NF}' $outsubdir/results-$imo.txt)
      OldError=$(awk '/^Error:/ {print $NF}' $outsubdir/results-$imo.txt)
      ErrorRate=$(perl -E "say sprintf(\"%.2f\", log($OldError/$Error)/log($OldMeshSize/$MeshSize))")
      echo -e "$Degree $MeshSize $NbCells $NbEdges $NbVertices $DimXDivDiv $Error $ErrorRate" >> $outsubdir/$errorsfile
    else
  	  echo -e "$Degree $MeshSize $NbCells $NbEdges $NbVertices $DimXDivDiv $Error --" >> $outsubdir/$errorsfile
    fi
done;
column -t < $outsubdir/$errorsfile | tee $outsubdir/$errorsfile
column -t < $outsubdir/$timesfile | tee $outsubdir/$timesfile

