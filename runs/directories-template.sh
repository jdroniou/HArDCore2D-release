# $origin is loaded in the runseries.sh script that calls directories.sh, and corresponds to the
# path of the folder from where this runseries.sh will be executed

# Source of all the project
root="/home/jdroniou/programmes/github/HArDCore2D"

# Location for the schemes' executable
executable=$root"/build/Schemes/$executable_name"

# Location of mesh files
meshdir=$root"/typ2_meshes"

# Location for all outputs. 
outdir=$origin"/outputs"

# Names for error file and latex file
errorsfile="data_rates.dat"
latexfile="rates.tex"

