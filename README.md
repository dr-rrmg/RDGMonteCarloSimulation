# Monte Carlo Simulation of Radon Collection in Different Gases

This repository contains the code and results for a Monte Carlo simulation of radon collection in different gases.

## Project Overview

This project aims to enhance our understanding of the collection efficiency of radon in various gases - a crucial aspect in research and industry. While experimental methods to derive this data exist, they are time-consuming and require access to various gases. This project presents an alternative - a Monte Carlo simulation that can predict a detector's response in different gases.

While there have been previous simulations of electrostatic collection, this project specifically models the widely-used DURRIDGE RAD7 device and accounts for chemical neutralisation processes, which are known to significantly impact collection.

## Electrostatic Collection Mechanism

The electrostatic collection process with RAD7 starts by sampling the gas of interest into the detector chamber. When 222Rn atoms enter the detector chamber, they distribute randomly as they are neutral and not responsive to the electric field. Post 222Rn alpha decay, a fraction of produced 218Po atoms becomes charged 218Po+ ions, which, if unneutralized, will continue their trajectory towards the detector due to the electric field. 

During its flight, 218Po+ is subject to collisions with carrier gas particles, which can cause deviations in its path. If the 218Po+ ion lands on the active surface of the detector, it is counted as a successful flight. If it lands elsewhere, it is counted as lost. It's important to note that the actual signal measured by the RAD7 is the alpha decay of the collected 218Po, which has a 50% chance of decaying into the detector's active surface.


## Modelling Electrostatic Collection Chamber E-Field

The electrostatic chamber and PIPS detector inside the RAD7 were modeled using ANSYS based on manufacturer dimensions. The chamber is a hemisphere on top of a cylinder, with a height of 6.1 cm, and a radius of 5.1 cm. The PIPS detector is placed along the cylinder axis on the flat side and has an active surface area of π(0.973 cm)².

The cylindrical symmetry of the electric field chamber was exploited to simplify the ANSYS electric field calculation from 3D to 2D by using cylindrical coordinates.

![Image description](images/RAD7Chamber.jpg)

Above are photos of the PIPS detector (left) and the internal chamber (right) in the RAD7.

ANSYS was used to simulate this geometry, which includes input parameters for material attributes and voltages. The walls and PIPS detector were modeled as metal, while the volume inside the chamber was modeled as gas. The chamber walls were set to a high voltage of 2500V, while the PIPS detector was set to 0V.

![Image description](path-to-image)

Above are ANSYS electric field solutions, on the left is a contour plot (in units of V/cm), and on the right is a vector direction plot with uniform arrow lengths. ANSYS generates an output file containing information about the magnitude and direction of the electric field for a mesh of nodes in the chamber. This electric field nodal map will be utilized in the simulation of 218Po+ drift.



## References

- Placeholder for references





