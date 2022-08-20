# optimal-sheeting-angle-control
Extremum Seeking Control (ESC) algorithms implementation for speed optimization of sail-assisted ships. The optimization is based on JAVAFOIL, an airfoil analysis software, that outputs the aerodynamic force. In this work, the thrust of the ship is maximized, and thus the output of interest from JAVAFOIL is the thurst coefficient $cT$.

The work was developed during an internship @ Chalmers University of Technology, in Sweden.

## Folders
- `\`: The root folder contains all the children folders as well as all the scripts.
- `lib\`: ESC implementation and interpolation scripts
- `data\`: Contains cT look-up tables, sheeting angle optimal values, and aparent wind angle (AWA) measurements from real experiments
- `docs\`: Contains the final report and presentation from this internship
- `plots\`: Contain plots from previous runs of the main scripts

The remaining folders are auxiliar to the JAVAFOIL calculations they are either well documented or its role in the software is easy to understand.

## Run
To run simulations, choose one of the following scripts:
- `esc`: Multivariate ESC implementation
- `esc_interp_`$x$`D`: Similar to `esc`, but with look-up table criterion allowing for faster computations. $x$ stands for number of dimensions of the schemes, i.e. number of sails in the simulation.
- `optimal_frequency_selection`: Script to obtain optimal frequency in the multivariate case. Runs a MILP program with the methods' constraints.