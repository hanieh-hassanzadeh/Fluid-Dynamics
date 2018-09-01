# Simulating Fluid Dynamics on Distributed Computers

The test-case selected here is to resolve the turbulent flow of a jet stream using [Large Eddy Simulation](https://en.wikipedia.org/wiki/Large_eddy_simulation) approach. In this approach the mean flow features are calculated by [Navier-Stokes](https://en.wikipedia.org/wiki/Navier–Stokes_equations) equations on a [grid](https://en.wikipedia.org/wiki/Grid_classification). The subgris scale features such as turbulence of the flow, however, are computed using parametrization equations. To solve the equations here, [Finite Difference](https://en.wikipedia.org/wiki/Finite_difference) disctitization is used. 

Although, LES approach helps to reduce the computational cost dramatically, it is still very costly to capture the turbulence of the flow. Therefore, two other approches were taken to decrease the computational cost; [implicit methods](https://en.wikipedia.org/wiki/Explicit_and_implicit_methods) and parallel processing to solve the equations.



* Disclosure: I have written this program in 2008. Both C++ and MPI versions used to compile this code are definitly dated. I will rewrite it with updated versions and commit it here. This version is for my own recored. However, some people may still find the cool idea I used to solve this problem interesting. 
