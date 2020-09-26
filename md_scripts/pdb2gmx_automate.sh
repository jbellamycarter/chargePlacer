#To be used with acp_gmx_command_list.sh script

gmx pdb2gmx -f protein_init.pdb -o protein.gro -p protein.top -ignh -ter -asp -glu -arg -lys -his -ff oplsaa -water none

echo "pdb2gmx has finished!"
