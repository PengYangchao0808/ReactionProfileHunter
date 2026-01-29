%NProcShared={nprocshared}
%Mem={mem}
%Chk={checkpoint_file}
#{method}/{basis} SCRF=(SMD,solvent={solvent}) EmpiricalDispersion={dispersion}
# IRC=(CalcFC, MaxPoints={max_points}, StepSize={step_size})

{charge} {multiplicity}
{ts_coordinates}

