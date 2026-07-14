# Finite Difference Schemes For Linear Advection Equation
This project explores finite difference schemes for solving the 1D Linear Advection Equation
Several numerical schemes were implemented and compared to analyze their accuracy, stability, and dispersion/dissipation characteristics.

## Implemented Schemes
The following schemes were coded and tested:
FTBS (Forward Time, Backward Space)
FTFS (Forward Time, Forward Space)
FTCS (Forward Time, Central Space)
Lax–Friedrichs Scheme
Lax–Wendroff Scheme
MacCormack Scheme
Crank–Nicolson Scheme

## Results & Comparisons
Each scheme was applied to the same initial condition.
Numerical solutions were plotted against the analytical solution.
Graphs highlight stability vs instability, numerical diffusion, and phase errors.

├── main.py              # Driver script
├── schemes/             # Implementations of each scheme
├── utils/               # Helper functions (plotting, error analysis)
├── results/             # Output graphs
└── README.md            # Project documentation

## Applications
Numerical PDE solving
Computational Fluid Dynamics (CFD)
Stability and accuracy analysis of numerical schemes

## Author
Developed by Kailash Bharat  
Computational Engineering, IIT Hyderabad
