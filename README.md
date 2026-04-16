# SemaGBA
A system dynamics model of the Semaglutide-responsive Gut-Brain Axis (SemaGBA). A model of how the brain and semaglutide regulate appetite and weight.

# Reproducing the Results
Run the following command from the `SemaGBA` folder to run all simulations and generate all available figures.

```bash
julia main.jl
```

This will run all the simulations and produce the figures from the paper and put them in the `figures` folder.

## Running Specific Simulations
To run specific simulations, use the following syntax:

```bash
julia main.jl 2 3
```

This example runs simulations 2 and 3, generating figures 2b and 3 respectively.

## Available Figures and Simulations
| Figure | Simulation | Description |
| --- | --- | --- |
| 1 | — | *Not generated in Julia* - Development workflow (created in BioRender) |
| 2a | — | *Not generated in Julia* - Causal loop diagram of reduced model (created in BioRender) |
| 2b | 'sim2b_Bar_plot.jl' | Bar plot of the model-predicted normalized change from baseline for metabolic variables following semaglutide treatment in type 2 diabetes (0.5 mg for 30 weeks) and obesity (2.4 mg for 68 weeks). |
| 3 | 'sim3_Overeating_plot.jl'| Model simulated progression from a healthy state to obesity with type 2 diabetes under chronic overeating with and without semaglutide intervention (2700 instead of 2,000 kcal/day) over 5 years |
| 4a | — | *Not generated in Julia* - Causal loop diagram of extended model (created in BioRender) |
| 4b | 'sim4b_Neuron_activity_plot.jl'|Model-predicted AgRP, POMC, and dopamine neural activity during semaglutide treatment in diabetes (0.5 mg semaglutide for 30 weeks) and obesity (2.4 mg semaglutide for 68 weeks) scenarios.|
| 5 | 'sim5_Diabetes_trajectories.jl'| Model simulated metabolic trajectories for diabetic patient treated with 0.5 mg semaglutide during 30 weeks (reduced and extended model comparison). |
| 6 | 'sim6_Obesity_trajectories.jl' |  Model simulated metabolic trajectories for obese patient treated with 2.4 mg semaglutide during 68 weeks (reduced and extended model comparison). |
| 7 | 'sim7_effect_functions_POMC.jl'| effect functions illustrating how insulin, leptin, GLP-1, AgRP and dopamine neural activity individually influence POMC neuron activity|
| 8 | 'sim8_contribution_to_neurons_plot.jl'| Contribution of GLP-1, leptin, insulin, dopamine neuron activity, and AgRP neuron activity to changes in AgRP, POMC, and dopamine neural activity over time in diabetes and obesity simulations.|

# Structure of this repository
## Needed Files for Simulations

### Core Model Files: src
- **Data.jl**: Contains the data points describing the effects of variables on each other, including:
  - Sigmoid function definition 
  - x- and y-data points describing the effects of variables on each other

- **Functions.jl**: Contains all utility functions required to run the simulations, including:
  - The main `run_simulation()` function
  - ODE callbacks for semaglutide and overeating interventions
  - Calculation functions for target variable values
  - Helper functions for interpolation and visualization
  
- **ODE_reduced.jl**: Defines the ODE system for the reduced model
- **ODE_extended.jl**: Defines the ODE system for the extended model (detailed neuronal dynamics)

### Input Files
- **Inputs_diabetes.jl**: Input parameters for diabetes simulations
- **Inputs_obesity.jl**: Input parameters for obesity simulations
- **Inputs_healthy.jl**: Input parameters for healthy baseline simulations

Each input file contains:
- Duration of the simulation (`t_days`)
- Pharmacokinetic parameters of semaglutide (bioavailability, distribution volume, absorption and elimination rates)
- Semaglutide injection doses and timepoints
- Reference, initial, and adjustment time values for all state variables

### Simulation Scripts
The simulations described in the table above generate the corresponding figures. Each script:
- Loads appropriate input files
- Runs one or more ODE solutions
- Computes target metrics or visualizations
- Saves the result as a PNG figure in the `figures` folder 

