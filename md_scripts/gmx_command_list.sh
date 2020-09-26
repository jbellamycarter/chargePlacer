# Gromacs Commands to run gas-phase MD Simulations with custom proton placement
# Jedd Bellamy-Carter, Apr 2020

# Prepare .gro and .top files, assumes file name is protein_init.pdb
# gmx pdb2gmx -f protein_init.pdb -o protein.gro -p protein.top -ignh -ter -asp -glu -arg -lys -his -ff oplsaa -water none

echo "$(tput setaf 2)Preparing input topology files$(tput sgr0)"
tput setaf 1
expect autocharger.exp

echo "$(tput setaf 2)Enlarging simulation box$(tput sgr0)"
# Make box big enough for unfolding and cut-offs
tput setaf 3
gmx editconf -f protein.gro -o protein_box.gro -c -box 900 900 900

echo "$(tput setaf 2)Minimising structure$(tput sgr0)"
tput setaf 4
# Perform energy minimisation of structure in-vacuo, requires min.mdp
gmx grompp -v -f min.mdp -c protein_box.gro -p protein.top -o min.tpr
gmx mdrun -v -deffnm min -c min.gro
echo Potential | gmx energy -f min.edr -o min_energy

gnuplot minim_plot.gnu

# Minimised structure is taken forward for all subsequent steps

# Equilibration and Production to be repeated 3 times with a different pseudorandom seed. i.e. from `gen-seed = -1`
#
# Perform equilibration runs for 50 ps

echo "$(tput setaf 2)Equilibrating structures$(tput sgr0)"
tput setaf 5
ROOT_DIR=$(pwd)
for DIR in RT A B C
do
    mkdir $ROOT_DIR/$DIR
    gmx grompp -v -f eq.mdp -c min.gro -r min.gro -p protein.top -o $DIR/eq.tpr
done

for DIR in RT A B C
do
    cd $ROOT_DIR/$DIR
    gmx mdrun -v -deffnm eq -s eq.tpr
    printf "Coulomb-(SR)\nPotential\nTemperature" | gmx energy -f eq.edr -o eq_energy
    gnuplot ../equil_plot.gnu
done

echo "$(tput setaf 2)Running MD simulations$(tput sgr0)"
tput setaf 6
cd $ROOT_DIR
gmx grompp -v -f prod_298K.mdp -c RT/eq.gro -t RT/eq.cpt -r RT/eq.gro -p protein.top -o RT/prod.tpr

for DIR in A B C
do
    gmx grompp -v -f prod_800K.mdp -c $DIR/eq.gro -t $DIR/eq.cpt -r $DIR/eq.gro -p protein.top -o $DIR/prod.tpr
done

for DIR in RT A B C
do
    gmx mdrun -v -deffnm $DIR/prod
done

echo "$(tput setaf 2)MD Simulations Completed!$(tput sgr0)"

# Analyse Production Run
for DIR in RT A B C
do
    cd $ROOT_DIR/$DIR
    echo Protein | gmx sasa -s eq.gro -f prod.xtc -o area.xvg
    echo Protein | gmx gyrate -s eq.gro -f prod.xtc -o gyrate.xvg
    echo Protein Protein | gmx rms -s eq.gro -f prod.xtc -o rmsd.xvg
    gnuplot $ROOT_DIR/../analysis_plot.gnu
done
# Create comparison figures
cd $ROOT_DIR
gnuplot $ROOT_DIR/analysis_all_plot.gnu

### Processing In PyMol
# 1. Open PyMol
# 2. cd to directory with MD runs
# 3. run the following commands: 
#       load eq.gro, 750K_run
#       load_traj prod.xtc, 750K_run
