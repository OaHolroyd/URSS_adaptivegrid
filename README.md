# URSS_adaptivegrid

![example](examples/wackywaves2.gif)

The repo houses the code/examples for my summer 2020 URSS project at the University of Warwick on the use of adaptive grids to reduce CPU time for solving the Benney equation IVP:

![equation1a](https://latex.codecogs.com/gif.latex?h_t%20-%20f%20&plus;%20q_x%20%3D%200)

![equation1b](https://latex.codecogs.com/gif.latex?q%20%3D%20%5Cfrac%7B2h%5E3%7D%7B3%7D%20-%20%5Cfrac%7Bh%5E3%7D%7B3%7D%5Cleft%28%202h_x%5Ccot%28%5Ctheta%29%20-%20%5Cfrac%7Bh_%7Bxxx%7D%7D%7BC%7D%5Cright%29%20&plus;%20R%20%5Cleft%28%20%5Cfrac%7B8h%5E6h_x%7D%7B15%7D%20-%20%5Cfrac%7B2h%5E4f%7D%7B3%7D%20%5Cright%29).

All the MATLAB code is in the [code](code) folder, which is further split into [solvers](code/solvers), containing the solvers and helper functions and [support](code/support), containing a variety of supporting functions including plotting scripts. All the functions and scripts are self-documented.

Example plots and gifs are all in the [examples](examples) folder. Most of these also include the `.mat` files required to reproduce them with the relevant plotting script.
