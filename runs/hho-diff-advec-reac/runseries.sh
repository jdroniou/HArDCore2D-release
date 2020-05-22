#!/bin/bash
#
# Execute hho file on series of meshes, and calculate outputs
#

###
# Executable file for the considered scheme
executable_name="hho-diff-advec-reac"

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
Executed without parameters: uses the data in the local data.sh file\n"
	exit;
fi;

# Create/clean output directory
if [ -d $outdir ]; then
	\rm -r $outdir
fi
	mkdir $outdir

###
# LOAD DATA
datafile=$(readlink -f data.sh);
echo -e "Using data file $datafile\n"
. $datafile

echo "degrees: (edge) k = $k, (cell) l = $l"
echo "boundary conditions bc = $bc"
echo -e "test case: solution = $tcsol, diffusion = $tcdiff, advection = $tcadvec, reaction = $tcreac"

###
# EXECUTE FOR EACH MESH
nbmesh=${#mesh[@]}
for i in `seq 1 $nbmesh`; 
do
  meshfile=${mesh[$i]}
  echo -e "\n*************************\nmesh $i out of $nbmesh: ${mesh[$i]}"
	# Execute code
  	$executable --mesh $meshfile --bc $bc --testcase $tcsol $tcdiff $tcadvec $tcreac --edgedegree $k  --celldegree $l --export_matrix $export_matrix --use_threads $use_threads 
	# Move outputs
	mv results.txt $outdir/results-$i.txt
	if [ -f T-"$plotfile".vtu ]; then
  	mv T-"$plotfile".vtu $outdir/mesh"$i"_T-"$plotfile".vtu
	fi
	if [ -f F-"$plotfile".vtu ]; then
  	mv F-"$plotfile".vtu $outdir/mesh"$i"_F-"$plotfile".vtu
  	fi
	if [ -f exact-$plotfile.vtu ]; then
  	mv exact-$plotfile.vtu $outdir/mesh"$i"_exact-$plotfile.vtu
  	fi
done;

# CREATE DATA FILE FOR LATEX
cd $outdir
echo -e "MeshSize L2Error H1Error EnergyError" > $errorsfile
for i in `seq 1 $nbmesh`; 
do	
	MeshSize=$(awk '/MeshSize:/ {print $NF}' results-$i.txt)
	L2Error=$(awk '/L2Error:/ {print $NF}' results-$i.txt)
	H1Error=$(awk '/H1Error:/ {print $NF}' results-$i.txt)
	EnergyError=$(awk '/EnergyError:/ {print $NF}' results-$i.txt)	
echo -e "$MeshSize $L2Error $H1Error $EnergyError" >> $errorsfile
done;

## CREATE AND COMPILE LATEX

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
        \\\addplot table[x=MeshSize,y=H1Error] {data_rates.dat};
        \\\addplot table[x=MeshSize,y=L2Error] {data_rates.dat};
        \\\logLogSlopeTriangle{0.90}{0.4}{0.1}{$kp2}{black};
        \\\logLogSlopeTriangle{0.90}{0.4}{0.1}{$kp1}{black};
        \\\legend{H1,L2};
      \\\end{loglogaxis}
    \\\end{tikzpicture}
  \\\caption{Error vs. \$h\$}
\\\end{figure}

Degrees: (edge) \$k=$k\$, (cell) \$l=$k\$

Boundary conditions: bc = $bc

Test case: solution = $tcsol, diffusion = $tcdiff, advection = $tcadvec, reaction = $tcreac \n\n" > $latexfile;

for i in `seq 1 $nbmesh`; 
do
	echo -e "mesh[$i] = \\\verb!${mesh[$i]}!\n\n" >> $latexfile;
done


echo -e "\\\end{document}" >> $latexfile;

# Tests data at the end of the file
for i in `seq 1 $nbmesh`; 
do
	echo -e "Test $i:\n" >> $latexfile
	cat results-$i.txt >> $latexfile
	echo -e "\n" >> $latexfile
	pdflatex $latexfile > /dev/null
done

##
# COMPUTATION OF CONVERGENCE RATES
#
echo -e "\n ----------- Data -----------"
echo "degrees: (edge) k = $k, (cell) l = $l"
echo "boundary conditions: bc = $bc"
echo "test case: solution = $tcsol, diffusion = $tcdiff, advection = $tcadvec, reaction = $tcreac"
cd ..
echo "data file: $datafile" >> $outdir/allrates.txt
echo "degrees: (edge) k = $k, (cell) l = $l" >> $outdir/allrates.txt
echo "boundary conditions: bc = $bc" >> $outdir/allrates.txt
echo "test case: solution = $tcsol, diffusion = $tcdiff, advection = $tcadvec, reaction = $tcreac" >> $outdir/allrates.txt
echo "Mesh 1: ${mesh[1]}" >> $outdir/allrates.txt
if [ ! -f ./compute_rates ]; then
  g++ compute_rates.cpp -o compute_rates
fi
./compute_rates | tee -a $outdir/allrates.txt
