# PARACHUTE++ 

PARACHUTE stands for: PARallel CHannel simUlaTor Evolution.

Realistic wireless channel models are an important part of any wireless network simulator. It brings a solid modelling of the physical world into the simulations and allows meaningful insights about techniques and networks performance and optimizations. 
With the evolution of wireless networks, the densification of its structure and the introduction of new wireless devices, the need of a performative channel generation library is adamant for research in wireless technologies.
In this context, this repository contains the implementation of continouos and discrete channel models and tools for channel modelling and analysis.
The reference material for this implementations is Mobile Radio Channels by Matthias Pätzold.
All the library is written in C++20, it is fully documented using doxygen and it is also based on C++ templates allowing to the user to chose the precision of numerical representation which can impact in storage or speed of calculations.

## Implemetation characteristics

- Fully C++ 20 compatible.
- Built in the OOP paradign.
- Optimized numerical algorithms.
- Fully documented with Doxygen.
- Uses STL containers for storage.
- GUI built in QT.
- Examples.
- Results tested against the theoretical models and the implemetation reference.

## Features

- Continuous Sum-of-Sinusoids models for time variant process generation.
- 6 different models for Sum-of-Sinusoids parametrization.
- Discrete random numbers generation with:
  - Rayleigh distribution.
  - Rice distribution.
  - lognormal distribution.
  - Suziki distribution.
  - gaussian distribution.
- Suzuki correlated models of Type I and II.
- CUDA-accelerated implementation of all models mentioned above.
