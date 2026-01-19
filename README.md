## SemaGBA
 A system dynamics model of the Semaglutide-responsive Gut-Brain Axis (SemaGBA). A model of how the brain and semaglutide regulate appetite and weight.

## Structure of this repository
# Input fils
- Input_diabetes.jl
- Input_obesity.jl
- Input_healthy.jl

These are the input files for diabetes, obesity and healthy simulations. Each file contains: 
- Duration of simulation (t_days)
- Pharmacokinetic parameters of semaglutide: bioavailability, distribution volume, absorption rate constant, and elimination rate constant. 
- Injection doses (injection_doses) and timepoints (release_times)
- Reference values (ref_), initial value (initial_) and adjustment time (_adj_time) of all variables. 

# Needed files for simulations
- Data.jl: Dontains the data points deeded to determine the effects of variables on each other. These data points are interpolated in the function file. 
- Functions.jl: Dontains all functions required to run the simulations, including the main run_simulation functio, ODE callbacks for semaglutide and overeating, and calculation of the target values. 
- ODE_simplified.jl: Defines the ODE system of the simplified model. 
- ODE_neural.jl: Defines the ODE system of the neural model.

# Simulation scripts for generating plots
- Bar_plot.jl: Generates bar plots showing normalized changes in variables relative to baseline for diabes and obesity simulations. 
- Overeating_plot.jl: Generates plots how variables change over time during chronic constant overeating, with and without semaglutide intervention. 
- Neural_activity_plot.jl: Generates plots of the neural activity (AgRP, POMC, and dopamine) over time during semaglutide treatment in diabetes and obesity simulations.  
- Diabetes_trajectories.jl: Generates metabolic trajectories for diabetes simulations. Both for the simplified and neural model. 
- Obesity_tjactories.jl: Generates metabolic trajectories for obesity simulation. Both for the simplified and neural model. 

