## Overview
This MATLAB simulation models the trajectory of a particle falling into a rotating black hole, described by the Kerr metric. The simulation accounts for gravitational effects, electric potential, and the conceptual aspect of Hawking radiation. It provides insights into the dynamics of particles near black holes.

## Features
- **Numerical Integration**: Uses the `ode45` function to solve the equations of motion for a particle influenced by gravitational and electric forces.
- **Hawking Radiation**: Conceptually incorporates the effects of Hawking radiation on the particle's energy.
- **Visualization**: Generates plots to visualize the radius, energy, angular momentum, and proper time of the particle throughout the simulation, as well as a 3D trajectory plot.

## Installation
To run this simulation, ensure you have MATLAB installed. No additional toolboxes are required.

## Usage
1. Open the `Plasma_Infused_RayBurst_Detonator_Railgun.m` file in MATLAB.
2. Run the script to execute the simulation.
3. View the output in the command window and the generated plots.

## Parameters
- **Black Hole Parameters**:
  - `M`: Mass of the black hole (arbitrary units).
  - `a`: Spin parameter (angular momentum per unit mass).

- **Particle Parameters**:
  - `m0`: Rest mass of the particle (arbitrary units).
  - `charge`: Charge of the particle (in electron charge units).
  - `potential`: Electric potential near the black hole (in volts).

- **Numerical Integration**:
  - `t_span`: Time span for integration (0 to 100 seconds).
  - `initial_conditions`: Initial conditions for the particle's motion.

## Results
The simulation outputs the final energy, angular momentum, and proper time experienced by the particle. It also provides energy values at different time points.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments
This simulation is based on concepts from general relativity and theoretical physics. For further reading, refer to literature on black holes, the Kerr metric, and Hawking radiation.
