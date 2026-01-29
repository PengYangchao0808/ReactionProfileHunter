%NProcShared={nprocshared}
%Mem={mem}
%Chk={checkpoint_file}
#{method}/{basis} SCRF=(SMD,solvent={solvent}) EmpiricalDispersion={dispersion}
# Opt=(QST2, CalcFC, NoEigenTest) Freq

{charge} {multiplicity}
{reactant_coords}

{charge} {multiplicity}
{product_coords}

