# Instructions for running simulation using the YLZ-poly potential

## Installation of LAMMPS with the required in-house changes

- Download the LAMMPS source code.
- Copy the custom source files (`YLZ-poly.cpp`, `YLZ-poly.h`, `add_ylz_pressure.cpp`, `add_ylz_pressure.h`, `compute_tension.cpp`, and `compute_tension.h`) into the `src/ASPHERE` directory.

- Compile and install LAMMPS (as outlined in LAMMPS documentation)
> Compiling with an MPI library for parallelization is highly recommended.

## Installation of the required Python libraries
- You need a python3 distribution with the following libraries
    -   numpy
    -   scipy
    -   ovito

## Running of the LAMMPS simulation
- Make sure that you have the `*.in` and `read_date.*` files in the same directory
> If you have already created a model in a `*.in` file, you don't need to read the `read_date.*` file.
- Execute the simulation using MPI: -in in.*`
  ```bash
  mpirun -np <number_of_cores> lmp_mpi -in in.flat_membrane
  ```
- The ouput files are:
    -   `log.lammps` = Standard LAMMPS log file
    -   `msd*.txt` = MSD (Mean Squared Displacement) output file
    -   `*.lammpstrj` = Dump output of particle trajectories
## Extracting Diffusion Coefficient from the Simulation
- Run Run the following Python script:
  ```bash
  python3 calculate_diffusion.py --file msd_results.txt
  ```

## Extracting fluctuation spectrum for the given simulation

- For this part, we use MATLAB.

- Ensure the `*.lammpstrj` file is placed in the same directory as your MATLAB scripts.
- Run the `helfrich_fitting_model.m` script.
