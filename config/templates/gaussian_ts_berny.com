%NProcShared={nprocshared}
%Mem={mem}
%Chk={checkpoint_file}
#{method}/{basis} SCRF=(SMD,solvent={solvent}) EmpiricalDispersion={dispersion}
# Opt=(TS, CalcFC, NoEigenTest) Freq

--
{charge} {multiplicity}
{coordinates}

--link1--
%NProcShared={nprocshared}
%Mem={mem}
%Chk={checkpoint_file}
#{method}/{basis} Geom=Check Guess=Read
# Freq=ReadIsotopes

--
