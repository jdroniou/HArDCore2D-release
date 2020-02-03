#!/bin/bash
#
# Execute hho file on series of meshes, and calculate outputs
#

###
# Executable file for the considered scheme
executable_name="hho-diffusion"

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
	echo -e "\nExecute tests using parameters in data.sh, creates and compile latex file, and calculate rates.\n
Executed without parameters: uses the data in the local data.sh file\n
Parameters that can be passed: data file | -k <value of k> | -l <value of l>\n"
	exit;
fi;

# Create/clean output directory
if [ -d $outdir ]; then
	\rm -r $outdir
fi
	mkdir $outdir

###
# LOAD DATA
# (LATER: Test that each required data exists (files, parameters...))
# Default data file
datafile=$(readlink -f data.sh);
# We go over the parameters. -k K or -l L will set k or l to these values. Any other
# parameter is assumed to be a datafile that overrides the previous one
until [ $# -eq 0 ]; do
  if [[ $1 == "-k" ]]; then
    ksave=$(echo $2 | sed s/'^[^0-9].*'//g)
    shift 2
  elif [[ $1 == "-l" ]]; then
    lsave=$(echo $2 | sed s/'^[^0-9].*'//g)
    shift 2
  else
    datafile=$(readlink -f $1)
    shift
  fi
done
echo -e "Using data file $datafile\n"
. $datafile

# if ksave or lsave have been ecountered we override whatever values $datafile might have contained
if [[ ! -z "$ksave" ]]; then
  k=$ksave;
fi
if [[ ! -z "$lsave" ]]; then
  l=$lsave;
fi

echo "degrees: (edge) k=$k, (cell) l=$l"
echo "boundary conditions bc=$bc"
echo -e "test case: solution=$tcsol, diffusion=$tcdiff\n"

###
# EXECUTE FOR EACH MESH
nbmesh=${#mesh[@]}
for i in `seq 1 $nbmesh`; 
do
  meshfile=$meshdir"/"${mesh[$i]}".typ2"
  echo -e "\n*************************\nmesh $i out of $nbmesh: ${mesh[$i]}.typ2"
	# Execute code
  $executable -m $meshfile -k $k -l $l -b $bc -c $tcsol $tcdiff --solver_type $solver_type --use_threads $use_threads --export_matrix $export_matrix
	r=$?
	if [ "$r" != "0" ]; then
		exit
	fi
	# Move outputs
	mv results.txt $outdir/results-$i.txt
	if [ -f T-solution.vtu ]; then
  	mv T-solution.vtu $outdir/mesh"$i"_T-solution.vtu
	fi
	if [ -f solution.vtu ]; then
  	mv solution.vtu $outdir/mesh"$i"_solution.vtu
  fi
	if [ -f exact-solution.vtu ]; then
  	mv exact-solution.vtu $outdir/mesh"$i"_exact-solution.vtu
  fi
	if [ -f SystemMatrix.mtx ]; then
  	mv SystemMatrix.mtx $outdir/mesh"$i"_SystemMatrix.mtx
  fi
done;

# CREATE DATA FILE FOR LATEX
cd $outdir
echo -e "meshsize L2error H1error EnergyError NbEdgeDOFs MeshReg MeshSkewness" > $errorsfile
for i in `seq 1 $nbmesh`; 
do
	meshsize=$(awk '/MeshSize:/ {print $NF}' results-$i.txt)
	L2error=$(awk '/L2error:/ {print $NF}' results-$i.txt)
	H1error=$(awk '/H1error:/ {print $NF}' results-$i.txt)
	EnergyError=$(awk '/EnergyError:/ {print $NF}' results-$i.txt)
  NbEdgeDOFs=$(awk '/NbEdgeDOFs:/ {print $NF}' results-$i.txt)
	MeshReg=$(awk '/MeshReg:/ {print $NF}' results-$i.txt)
  MeshSkewness=$(awk '/MeshSkewness:/ {print $NF}' results-$i.txt)
  echo -e "$meshsize $L2error $H1error $EnergyError $NbEdgeDOFs $MeshReg $MeshSkewness" >> $errorsfile
done;

## CREATE AND COMPILE LATEX

# Look for minimal meshsize, that we reduce a bit
a=$(cat $errorsfile | awk ' END{print} ' | cut -d ' ' -f 1);
xmin=$(echo ${a/[eE]/*10^}*0.6 | bc -l);

kp2=$(($k+2));
kp1=$(($k+1));

echo -e "\\\documentclass{article}

\\\usepackage{amsfonts,latexsym,graphicx}

\\\setlength{\\\textwidth}{16cm}
\\\setlength{\\\textheight}{23cm}
\\\setlength{\\\oddsidemargin}{0cm}
\\\setlength{\\\evensidemargin}{0cm}
\\\setlength{\\\topmargin}{-1cm}
\\\parindent=0pt

\\\usepackage{pgfplots,pgfplotstable}
\\\usetikzlibrary{calc,external}

\\\newcommand{\\\logLogSlopeTriangle}[5]
{
    \pgfplotsextra
    {
        \\\pgfkeysgetvalue{/pgfplots/xmin}{\\\xmin}
        \\\pgfkeysgetvalue{/pgfplots/xmax}{\\\xmax}
        \\\pgfkeysgetvalue{/pgfplots/ymin}{\\\ymin}
        \\\pgfkeysgetvalue{/pgfplots/ymax}{\\\ymax}

        % Calculate auxilliary quantities, in relative sense.
        \\\pgfmathsetmacro{\\\xArel}{#1}
        \\\pgfmathsetmacro{\\\yArel}{#3}
        \\\pgfmathsetmacro{\\\xBrel}{#1-#2}
        \\\pgfmathsetmacro{\\\yBrel}{\\\yArel}
        \\\pgfmathsetmacro{\\\xCrel}{\\\xArel}

        \\\pgfmathsetmacro{\\\lnxB}{\\\xmin*(1-(#1-#2))+\xmax*(#1-#2)} % in [xmin,xmax].
        \\\pgfmathsetmacro{\\\lnxA}{\\\xmin*(1-#1)+\xmax*#1} % in [xmin,xmax].
        \\\pgfmathsetmacro{\\\lnyA}{\\\ymin*(1-#3)+\ymax*#3} % in [ymin,ymax].
        \\\pgfmathsetmacro{\\\lnyC}{\\\lnyA+#4*(\\\lnxA-\\\lnxB)}
        \\\pgfmathsetmacro{\\\yCrel}{\\\lnyC-\\\ymin)/(\\\ymax-\\\ymin)}

        % Define coordinates for \draw. MIND THE 'rel axis cs' as opposed to the 'axis cs'.
        \\\coordinate (A) at (rel axis cs:\\\xArel,\\\yArel);
        \\\coordinate (B) at (rel axis cs:\\\xBrel,\\\yBrel);
        \\\coordinate (C) at (rel axis cs:\\\xCrel,\\\yCrel);

        % Draw slope triangle.
        \\\draw[#5]   (A)-- node[pos=0.5,anchor=north] {\\\scriptsize{1}}
                    (B)-- 
                    (C)-- node[pos=0.,anchor=west] {\\\scriptsize{#4}} %% node[pos=0.5,anchor=west] {#4}
                    cycle;
    }
}

\\\begin{document}

\\\begin{figure}\\\centering
  \\\ref{legend}
  \\\vspace{0.5cm}\\\\\\
    \\\begin{tikzpicture}
      \\\begin{loglogaxis}[ legend columns=-1, legend to name=legend ]
        \\\addplot table[x=meshsize,y=EnergyError] {data_rates.dat};
        \\\addplot table[x=meshsize,y=H1error] {data_rates.dat};
        \\\addplot table[x=meshsize,y=L2error] {data_rates.dat};
        \\\logLogSlopeTriangle{0.90}{0.4}{0.1}{$kp2}{black};
        \\\logLogSlopeTriangle{0.90}{0.4}{0.1}{$kp1}{black};
        \\\legend{Energy,H1,L2};
      \\\end{loglogaxis}
    \\\end{tikzpicture}
  \\\caption{Error vs. \$h\$}
\\\end{figure}

Degrees: (edge) \$k=$k\$, (cell) \$l=$l\$

Boundary conditions: bc=$bc

Test case: tcsol=$tcsol, tcdiff=$tcdiff \n\n" > $latexfile;

for i in `seq 1 $nbmesh`; 
do
	echo -e "mesh[$i]=\\\verb!${mesh[$i]}!\n\n" >> $latexfile;
done


echo -e "\\\end{document}" >> $latexfile;


# Tests data at the end of the file
for i in `seq 1 $nbmesh`; 
do
	echo -e "Test $i:\n" >> $latexfile
	cat results-$i.txt >> $latexfile
	pdflatex $latexfile > /dev/null
done

##
# COMPUTATION OF CONVERGENCE RATES
#
echo -e "\n ----------- Data -----------"
echo "degrees: (edge) k=$k, (cell) l=$l"
echo "boundary conditions bc=$bc"
echo "test case: solution=$tcsol, diffusion=$tcdiff"
cd ..
echo "data file: $datafile" > $outdir/allrates.txt
echo "k=$k, l=$l" >> $outdir/allrates.txt
echo "Mesh 1: ${mesh[1]}" >> $outdir/allrates.txt
make compute_rates_run | tee -a $outdir/allrates.txt


