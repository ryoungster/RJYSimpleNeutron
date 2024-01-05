SimpleNeutron

INTRO:
Implementing Monte Carlo and Finite-Difference neutron transport/diffusion numerical models as a header based library in C. 
Based on course-work for NSE 451.

USAGE:
Presently only Gauss-Seidel Finite Difference is implemented for steady state, non-multiplying slabs (Plus other approximations in diffusion). To use:
    Include "Neutronic.h"
    Format your region data into the RegionDesc struct:
        Region Count
        Arrays of Source Strength, Absorbtion Macroscopic Cross-Section, Diffusion Coefficient, Region Length per Region
        Array of Cells per Region
    Run "GaussSeidel(RegionDesc, Epsilon)" with your chosen covergence criteria
    Returned data will be a GSOutput struct with with members:
        Cell count (Array Lengths are 1+ this value due to edges)
        Array of X-Values
        Array of Phi-Values
    In the examples this data is graphed using DISLIN (www.dislin.de) a free graphing library.

Given compilation script for examples expects DISLIN libraries and windows.

PLANS:
Refactor code slightly
Re-evaluate choice to put everything in one header
Power iteration/Monte Carlo k-calc
Standardise Input formating into something nice (look into dataset based internal calculation)
Time Dependent Calculations

POSSIBLE PLANS:
Error handling (As in malloc fails, etc.)
Multithreading 
