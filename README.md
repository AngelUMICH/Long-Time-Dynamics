# Long-Time-Dynamics (OPLD)
PhD Research Repository of Angel Chavira at the University of Michigan NERS department.

## Background

OPLD is a code which was designed to assist in the study of cascade overlap and defect annealing on damage and time scales up to 1 dpa and a few seconds,
respectively. Irradiation cascades are modeled by two methods: random frenkel pair insertion analogous to electron irradiation and actual ion damage events
whose dynamics are driven by traditional MD simulation. With atomic precision defects created by the damage events can be modeled by interatomic forces and 
subsequently their transition states can be found by the dimer search.

## Building Code

### Requirements

The LAMMPS **static** library must be built with the support of OpenMPI and must include any necessary packages needed to run MD such as the MANYBODY pakage (enables support for eam force styles). After the static library is built the path to the library must be included in the Makefile as well as the LAMMPS code source.

### Next Steps

Traditional make is used to create the binaries and the executables. Avoid using the multi-threaded make approach as this way of building is not consistent and might fail.

### How to Use OPLD

OPLD can be initiated two different ways: 1) new simulation or 2) restarted simulation. The former method of using OPLD requires only two input files: 'par_init.dat' and 'in.lammps'. Whereas, the latter method of using the code requires an additional file called 'RESTART.OUT'. As a note, the files must be included in the working directory with the correct name in order for the program to find them.

## OPLD Input file

An input file named 'par_init.dat' must be availiable for the code to read in. Below is a sample of what the input file must look like. The ordering of 
parameters is important. The first column contains the values of the parameters, second columns are the names of the variables used in the code. In the third
column within the brackets are the units of each parameters. To the right of the brackets is the data type of the parameters. Finally, the last column is a
short description of each of the parameters.

    1           # REST          [ --  ] I  Using Restart File? 1:No 2:Yes
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

## LAMMPS Input File
The file named 'in.lammps' is one of the two files that must always be present when running any kind of simulation. This file is used to set various features within the LAMMPS ecosystem which is needed to calculate forces and energies of an entire system or a subset atomic volume. This that this file must always set is potential type and coefficients. Whenever, the program is not being used to restart a simulation (i.e not in restart mode) the LAMMPS input file must be used to specify the initial atomic coordinates. Through LAMMPS this can be done in various manners such as reading a LAMMPS data or restart file, or an XYZ file. Keep in mind that a LAMMPS restart file is not the same as the 'RESTART.OUT' which is an output of the OPLD software.

The last thing that this input file must contain is the mass of all active species which is necessary for running a basic MD simulation. All other commands such as computes, groups, fixes, etc. should not be included into the file. The OPLD program all ready handles all necessary comptues and fixes. Including the other LAMMPS commands which are not sanction or stated here might cause the program to misbehave or crash.

Below is a sample 'in.lammps' input file. In this example, the user has opted out from reading in a system using a restart file or a data file and instead has created the system from LAMMPS commands. The potential type and coefficients were also set. Finally, the user chose to change the timestep. Including this command will only have an effect when simulating cascades as all other parts of this code do not run dynamics.

    lattice bcc 2.8553
    region box block 0 30 0 30 0 30
    create_box 1 box
    create_atoms 1 box
    mass 1 56.0
    
    pair_style eam/fs
    pair_coeff * * FeP_mm.eam.fs Fe
    neighbor 0 bin
    
    timestep 0.002

When running this program in restart mode it is best to leave the LAMMPS input file the same as when the restart file was written. In other words, use the same 'in.lammps' file when using a 'RESTART.OUT' file. A restart file does not include information about the initial configuration of the system.

## Using a Restart File

A simulation can be restarted from specific points by use of a restart file. These restart files are created automatically and stored in the 'RESTART' directory. Restarting a simulation can be convinient because it can act as a checkpoint after a failed run. The information that is stored in the restart file are values and parameters which are not already stored in the 'par_init.dat' file. Only information that pertains to that specific simulation is stored such as migration energies for defects, statistical values pertaining to damage events, and atom positions.

In order to use a restart file the correct option in the 'par_init.dat' file must be set. Also, the restart file must be placed in the working directory and renamed to **RESTART.OUT**. Using a restart file should not change any already created data in the working directory.