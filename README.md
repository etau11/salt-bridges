# salt-bridges
a tcl script calculating salt bridges between two groups for each residue
```
vmd -dispdev text -e salt-bridges.tcl
```
To calculate the probability of residues in group A and group B potentially engaging in salt-bridge interactions, simply input the ranges of group A and group B (modify the "resid" in the script). The script will ouput the results and save in sb_probabilities.dat.
