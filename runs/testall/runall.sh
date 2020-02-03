#!/bin/bash
#
# Execute "runseries" for a all polynomial degrees up to kmax and lmax
#

dirall=$(pwd);
outdir=$dirall"/outputs/";


# Check parameter
if [[ $# == 0 ]]; then
  exit;
fi;
scheme=$1;
# Load data
datafile=$(readlink -f data.sh);
if [ $# == 2 ]; then
  datafile=$(readlink -f $2);
fi
echo -e "Using data file $datafile\n"
. $datafile

# Create/clean output directory
if [ -d $outdir ]; then
	\rm -r $outdir
fi
mkdir $outdir


# Go the scheme directory
cd ../$scheme;
pwd;

# Save data file in scheme's directory
cp data.sh save-data.sh;

# To determine the bounds of l: boundsl='all' if from k-1 to k+1, otherwise just l=k
#boundsl="all";
boundsl="k";

# Run the scheme for all values of k,l
for (( k=0; k<=$kmax; k++))
do
  # lower bounds and upper bounds for l (from k-1 to k+1)
  if [ $boundsl == "all" ]; then
    lowerl=$(($k-1));
    upperl=$(($k+1));
  else
    lowerl=$k;
    upperl=$k;
  fi

  for (( l=$lowerl; l<= $upperl; l++))
  do
    if [[ $l -le $lmax ]];  then
      echo "k=$k, l=$l";
      cat $datafile | sed s/'^kmax=.'/'k='$k/ | sed s/'^lmax=.'/'l='$l/ > data.sh;
      cp data.sh outputs/;
      ./runseries.sh > alloutput.k$k.l$l.txt;
      mv alloutput.k$k.l$l.txt outputs/;
      ./COMPUTE_RATES > rates.k$k.l$l.txt;
      cp -r outputs/ $outdir/k$k.l$l/;
      mv rates.k$k.l$l.txt $outdir;
  sleep 1;
    fi
  done
done

# Recover data.sh originally in the scheme directory
mv save-data.sh data.sh;

# Concatenate all rates
cd $outdir;
echo $scheme >> allrates.txt;
for (( k=0; k<=$kmax; k++))
do
  echo -e "\n\n########## k=$k ###########" >> allrates.txt;
  # lower bounds and upper bounds for l (from k-1 to k+1)
  if [ $boundsl == "all" ]; then
    lowerl=$(($k-1));
    upperl=$(($k+1));
  else
    lowerl=$k;
    upperl=$k;
  fi

  for (( l=$lowerl; l<= $upperl; l++))
  do
    if [[ $l -le $lmax ]];  then
      echo -e "\n ****** k=$k, l=$l ******" >> allrates.txt;
      cat rates.k$k.l$l.txt >> allrates.txt;
      echo -e "\n" >> allrates.txt;
    fi
done
done


