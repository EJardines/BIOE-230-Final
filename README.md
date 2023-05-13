# Final Project for PHYS-230/BIOE-230 by Eira Jardines and Jesselynn LaBelle

## Contents
This repository contains files related to the final project for PHYS-230/BIOE-230. To run the simulation, clone this repository, open the `Simulation_v07.m` file, and run the script. Some special packages may need to be downloaded in order to run the script, depending on what is already included on the user's computer.

## Simulation
The simulation is written in Matlab within the file entitled `Simulation_v07.m`. This simulation was built upon a peer's code to model Contact Inhibition of Locomotion.

## Important Edits
* The cell collision was removed from the code (by preventing a collision check to occur)
* The number of cells were reduced to only 10
* Replaced periodic boundaries with solid boundaries
* Simulated cells on the surface of half a hemisphere instead of a rectangular region (to simulate cells on a petri dish)
* Cells to start at the bottom of half the hemisphere
* Removed noise in the initial y-axis position
* Attractive force was increased (to simulate a reduction in repulsion)

## Post-Presentation Changes
After the class presentation, an attempt to create a solid circular boundary for the cells to bounce off of lead to many issues. An adjustment to the location of the hemisphere had to be made to allow for easier calculations. The collision between the cell and hemisphere was done by calculating the particle's radius at each position using the Pythagorean Theorem. If the radius was greater than or equal to the radius of the petri dish, the code altered the rest of the particle's position to simulate a collision. The bottom half of the hemisphere was also treated as a solid wall and flipped the particle's position if a collision was detected. There are still issues with the collision detection at the top of the hemisphere. 


# ReadMe is adapted from this previous project
# Final Project for PHYS-230/BIOE-230

## Contents
This repository contains files related to the final project for PHYS-230/BIOE-230. To run the simulation, clone this repository, open the `Simulation_v07.m` file, and run the script. Some special packages may need to be downloaded in order to run the script, depending on what is already included on the user's computer.

## Simulation
The simulation is written in Matlab within the file entitled `Simulation_v07.m`. Two coupled Langevin equations are used, one for two-dimensional translational motion and one for one-dimensional rotational motion. Periodic boundary conditions are enforced.

Currently, there are some redundant calculations within the time loop, but those can be easily fixed by anyone interested in doing so. The code is commented and should therefore be somewhat easy to follow for those who are already familiar with Euler methods and the Langevin equation.

## Experimental Data
Experimental data are included within the two Excel files named `ExperimentalData.xlsx` and `ExperimentalData_02.xlsx`. From this data, rotational and translational coefficients are extracted and plotted using the `ExperimentalDataAnalysis.m`, which calls the `ImportExcelFile.m` function, which was automatically constructed by Matlab to pull data from Excel sheets whose format is identical to that of `ExperimentalData.xlsx` and `ExperimentalData_02.xlsx`.

## Post-Presentation Changes
Following the class presentation, various attempts were made to produce simulation results that exhibited behavior similar to that shown in videos of zebrafish cells, such as collective motion. However, the suggested changes to the code were not enough to produce desired results. As such, rather than invest further time into rewriting the code, the decision was made to submit what was already complete, though we understand the results are less than desirable.
