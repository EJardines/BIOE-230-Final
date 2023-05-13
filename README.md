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
