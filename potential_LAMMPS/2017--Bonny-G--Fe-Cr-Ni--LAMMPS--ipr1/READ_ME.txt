The potentials in LAMMPS are read by the following three commands:

pair_style hybrid/overlay eam/alloy eam/fs
pair_coeff * * eam/alloy FeCrNi_d.eam.alloy Fe Cr Ni
pair_coeff * * eam/fs FeCrNi_s.eam.fs Fe Cr Ni