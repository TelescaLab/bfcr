# Bayesian Functional Covariance Regression

Please see README.ipynb.


## Running Simulations

The simulations conducted in this paper can be run using the following steps:

   1. Open a terminal window and naviage to the `/Simulations` folder
   2. To run the N = 100 simulation, type `sh run_sim_100.sh` into the command line and press enter
   3. To run the N = 400 simulation, type `sh run_sim_400.sh` into the command line and press enter

Note: To change the number of cpus used in the simulation, open up the `Run_simulations_x00.R` file and change `ncpu <- min(4, availableCores())`. To change the number of MCMC iterations (50,000) or the burn-in parameter (25,000 iterations), change `run_mcmc(..., 50000, 25000, 5)`.
