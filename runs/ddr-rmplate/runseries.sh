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
executable_name="ddr-rmplate"
#make $executable_name

origin=$(pwd)
if [ ! -f ../directories.sh ]; then
  echo "directories.sh does not exist. Please read the README.txt in the parent folder."
  exit
fi
. ../directories.sh

for t in $ts; do

  echo "Degree                 : $k"
  echo "Mesh family            : $mesh_family"
  echo "Stabilization parameter: $stab_par"
  echo "t/E/nu                 : $t/$E/$nu"
  echo "Test case: solution    : $solution"
  echo "Output directory       : $outdir"

  outsubdir=${outdir}/${mesh_family}_k${k}_t${t}
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
      if($executable -m $meshfile.typ2 -k $k -s $solution -b $BC -x $stab_par -t $t -E $E -n $nu -e $export_matrix); then
        # Move outputs
        mv results.txt $outsubdir/results-$i.txt
        if [ -f displacement.vtu ]; then
        	mv displacement.vtu $outsubdir/mesh"$i"_displacement.vtu
      	fi
      	if [ -f exact-displacement.vtu ]; then
        	mv exact-displacement.vtu $outsubdir/mesh"$i"_exact-displacement.vtu
        fi
      	if [ $export_matrix != "false" ]; then
        	mv A_rm.mtx $outsubdir/A"$i".mtx
        	mv b_rm.mtx $outsubdir/b"$i".mtx
        fi
      fi
  done

  # CREATE DATA FILE FOR LATEX
  echo -e "Deg MeshSize NbCells NbEdges DimEXCurl DimXGRad EnError Rate DisCentError Rate" > $outsubdir/$errorsfile
  echo -e "TwallDDRCore TprocDDRCore TwallModel TprocModel TwallSolve TprocSolve" > $outsubdir/$timesfile
  for i in `seq 1 $nbmesh`; 
  do
      Degree=$(awk '/Degree:/ {print $NF}' $outsubdir/results-$i.txt)
      MeshSize=$(awk '/MeshSize:/ {print $NF}' $outsubdir/results-$i.txt)
      NbCells=$(awk '/NbCells:/ {print $NF}' $outsubdir/results-$i.txt)
      NbEdges=$(awk '/NbEdges:/ {print $NF}' $outsubdir/results-$i.txt)
      DimEXCurl=$(awk '/DimEXCurl:/ {print $NF}' $outsubdir/results-$i.txt)
      DimXGrad=$(awk '/DimXGrad:/ {print $NF}' $outsubdir/results-$i.txt)
      EnError=$(awk '/EnError:/ {print $NF}' $outsubdir/results-$i.txt)
      DisCentError=$(awk '/DisCentError:/ {print $NF}' $outsubdir/results-$i.txt)
      TwallDDRCore=$(awk '/TwallDDRCore:/ {print $NF}' $outsubdir/results-$i.txt)
      TprocDDRCore=$(awk '/TprocDDRCore:/ {print $NF}' $outsubdir/results-$i.txt)
      TwallModel=$(awk '/TwallModel:/ {print $NF}' $outsubdir/results-$i.txt)
      TprocModel=$(awk '/TprocModel:/ {print $NF}' $outsubdir/results-$i.txt)
      TwallSolve=$(awk '/TwallSolve:/ {print $NF}' $outsubdir/results-$i.txt)
      TprocSolve=$(awk '/TprocSolve:/ {print $NF}' $outsubdir/results-$i.txt)
      echo -e "$TwallDDRCore $TprocDDRCore $TwallModel $TprocModel $TwallSolve $TprocSolve" >> $outsubdir/$timesfile
      if(($i > 1)); then
        imo=$(perl -E "say $i - 1")
        OldMeshSize=$(awk '/MeshSize:/ {print $NF}' $outsubdir/results-$imo.txt)
        OldEnError=$(awk '/EnError:/ {print $NF}' $outsubdir/results-$imo.txt)
        EnErrorRate=$(perl -E "say sprintf(\"%.2f\", log($OldEnError/$EnError)/log($OldMeshSize/$MeshSize))")
        OldDisCentError=$(awk '/DisCentError:/ {print $NF}' $outsubdir/results-$imo.txt)
        DisCentErrorRate=$(perl -E "say sprintf(\"%.2f\", log($OldDisCentError/$DisCentError)/log($OldMeshSize/$MeshSize))")
        echo -e "$Degree $MeshSize $NbCells $NbEdges $DimEXCurl $DimXGrad $EnError $EnErrorRate $DisCentError $DisCentErrorRate" >> $outsubdir/$errorsfile
      else
    	  echo -e "$Degree $MeshSize $NbCells $NbEdges $DimEXCurl $DimXGrad $EnError -- $DisCentError --" >> $outsubdir/$errorsfile
      fi
  done;
  column -t < $outsubdir/$errorsfile | tee $outsubdir/$errorsfile
  column -t < $outsubdir/$timesfile | tee $outsubdir/$timesfile

done;  # end for t in $ts
