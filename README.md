# SemaGBA
A system dynamics model of the Semaglutide-responsive Gut-Brain Axis (SemaGBA). A model of how the brain and semaglutide regulate appetite and weight.

# Reproducing the Results
Run the following command from the `SemaGBA` folder to run all experiments.

```bash
julia main.jl
```

This will run all the simulations and produce the figures from the paper and put them in the `figures` folder.

If you want to run specific experiments, use the following syntax:

```bash
julia main.jl 1 2
```

to run the experiments producing figures 1 and 2. 

| Experiment number | Description |
| --- | --- |
| 1 | beschrijving |


# Structure of this repository
## Input fils
- **Input_diabetes.jl**
- Input_obesity.jl
- Input_healthy.jl

These are the input files for diabetes, obesity and healthy simulations. Each file contains: 
- Duration of the simulation (t_days)
- Pharmacokinetic parameters of semaglutide: bioavailability, distribution volume, absorption rate constant, and elimination rate constant. 
- Semaglutide's injection doses (injection_doses) and timepoints (release_times)
- Reference values (ref_), initial value (initial_) and adjustment time (_adj_time) of all variables. 

## Needed files for simulations
- **Data.jl**: Contains the data points describing the effects of variables on each other. 
- **Functions.jl**: Contains all functions required to run the simulations, including:
    - The main run_simulation function
    - ODE callbacks for semaglutide and overeating
    - Calculation of the target values
    - Functions to generate the bar plot 
- **ODE_simplified.jl**: Defines the ODE system of the simplified model. 
- **ODE_neural.jl**: Defines the ODE system of the neural model.

## Simulation scripts for generating plots
- **Bar_plot.jl**: Generates bar plots showing normalized changes in variables relative to baseline for diabetes and obesity simulations. 
- **Overeating_plot.jl**: Generates plots how variables change over time during chronic constant overeating, with and without semaglutide intervention. 
- **Neural_activity_plot.jl**: Generates plots of the neural activity (AgRP, POMC, and dopamine) over time during semaglutide treatment in diabetes and obesity simulations.  
- **Diabetes_trajectories.jl**: Generates metabolic trajectories for diabetes simulations. Both for the simplified and neural model. 
- **Obesity_tjactories.jl**: Generates metabolic trajectories for obesity simulation. Both for the simplified and neural model. 

