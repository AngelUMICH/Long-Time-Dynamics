# Long-Time-Dynamics
PhD Research Repository of Angel Chavira at the University of Michigan NERS department.

## Background

OPLD is a code which was designed to assist in the study of cascade overlap and defect annealing on damage and time scales up to 1 dpa and a few seconds,
respectively. Irradiation cascades are modeled by two methods: random frenkel pair insertion analogous to electron irradiation and actual ion damage events
whose dynamics are driven by traditional MD simulation. With atomic precision defects created by the damage events can be modeled by interatomic forces and 
subsequently their transition states can be found by the dimer search.

## Sample Input file

An input file named 'par_init.dat' must be availiable for the code to read in. Below is a sample of what the input file must look like. The ordering of 
parameters is important. The first column contains the values of the parameters, second columns are the names of the variables used in the code. In the third
column within the brackets are the units of each parameters. To the right of the brackets is the data type of the parameters. Finally, the last column is a
short description of each of the parameters.

    2           # TJOB    	    [ --  ] I  Job Type	: 1:Relaxation; 2:OPLD
    2           # TDAM          [ --  ] I  Irradition mode : 1:Cascade; 2:Electron
    1           # TPBC          [ --  ] I  Boundary Condition: 1: PBC 2: Absorbing Boundary Condition
    3           # NTSTEP        [ --  ] I  Total Number of Step
    1           # NFREQ         [ --  ] I  Output Frequency For Ovito
    1500        # NDEMAX        [ --  ] I  Maximum Number of Defect
    7000        # NAVMAX        [ --  ] I  Maximum Number of Atoms in the AVs
    12          # NDIMER        [ --  ] I  Number of dimer searches per iteration 
    2.8553      # RLAT          [ A   ] R  Lattice Constant
    30.0        # RBOX          [ a   ] R  Side of the Simulation Box
    5.3         # PCUTR         [ A   ] R  Cutoff for the pair    term of potential
    4.21        # DCUTR         [ A   ] R  Cutoff for the density term of potential
    0.20        # DECUT         [ a   ] R  Cutoff for Find Defects Wigner-Ziets
    3           # NDRMAX        [ --  ] I  Maximum Number of Rotation in Each Dimer Step
    0.001       # RDIMER        [ a   ] R  The Distance Between the Dimer
    0.100       # FNRMAX        [ --  ] R  The Maximum Rotation Force
    0.001       # FNRMIN        [ --  ] R  The Minimum Rotation Force
    0.100       # DSTPEF        [ a   ] R  The Step Size in Dimer Translation using CG
    0.070       # DSTEPM        [ a   ] R  The Maximum Step Size in Dimer Translation
    0.800       # EMSKIP        [ eV  ] R  Maximum Acceptable Energy Barrier
    300.000     # TEMP          [ k   ] R  Temperature
    300.0       # PKAE          [ ev  ] R  PKA Energy (Ion Mode only)
    1.20        # DRATE         [dpa/s] R  Irradiation Dose Rate
    10          # NFREKL        [ dpa ] I  Frenkel Pairs Per Damage Sequence (Electron Mode Only)
    0.000E+00   # TIME          [ s   ] R  Time
