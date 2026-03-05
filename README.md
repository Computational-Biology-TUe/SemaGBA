# SemaGBA
A system dynamics model of the Semaglutide-responsive Gut-Brain Axis (SemaGBA). A model of how the brain and semaglutide regulate appetite and weight.

# Reproducing the Results
Run the following command from the `SemaGBA` folder to run all simulations.

```bash
julia main.jl
```

This will run all the simulations and produce the figures from the paper and put them in the `figures` folder.

If you want to run specific simulations, use the following syntax:

```bash
julia main.jl 2 3
```

to run the simulations producing figures 2 and 3. 

| Simulation number | Description |
| --- | --- |
| 2 | Figure 2; Bar plot of the normalized change from baseline for metabolic variables following semaglutide treatment in type 2 diabetes (0.5 mg for 30 weeks) and obesity (2.4 mg for 68 weeks). |
| 3 | Figure 3; Simulation of the progression from a healthy state to obesity with type 2 diabetes under chronic overeating (3000 instead of 2000 kcal/day) during 3 years. Comparison between untreated progression and treatment with 1.0 mg semaglutide treatment initiated at day 210. |
| 4 | Figure 4b; AgRP, POMC, and dopamine neuron activity in two model scenarios; diabetes (0.5 mg semaglutide for 30 weeks) and obesity (2.4 mg semaglutide for 68 weeks).|
| 5 | Figure 5; Metabolic trajectories of a diabetic patient treated with 0.5 mg semaglutide during 30 weeks, including the reduced and extended simulation. |
| 6 | Figure 6; Metabolic trajectories of an obese patient treated with 2.4 mg semaglutide during 68 weeks, including the reduced and extended simulation.|



# Structure of this repository
## Input fils
- **Input_diabetes.jl**
- **Input_obesity.jl**
- **Input_healthy.jl**

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
- **ODE_reduced.jl**: Defines the ODE system of the reduced model. 
- **ODE_extended.jl**: Defines the ODE system of the extended model.

## Simulation scripts for generating plots
- **Bar_plot.jl**: Generates bar plots showing normalized changes in variables relative to baseline for diabetes and obesity simulations. 
- **Overeating_plot.jl**: Generates plots how variables change over time during chronic constant overeating, with and without semaglutide intervention. 
- **Neural_activity_plot.jl**: Generates plots of the neural activity (AgRP, POMC, and dopamine) over time during semaglutide treatment in diabetes and obesity simulations.  
- **Diabetes_trajectories.jl**: Generates metabolic trajectories for diabetes simulations. Both for the reduced and extended model. 
- **Obesity_tjactories.jl**: Generates metabolic trajectories for obesity simulation. Both for the reduced and extended model. 

