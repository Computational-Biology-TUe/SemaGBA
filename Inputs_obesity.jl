""" This file contains all inputs to simulate an obese patient"""

# Define the duration of the simulation
t_days = 61*7; #days

# Define the pharmacokinetic parameters for the semaglutide injections
bioavailability = 0.89; #dimentionless
distribution_volume = 125.0; #dL

# Define elimination and absorption rate of semaglutide
absorption_rate = 0.0286 *24 #/day
elimination_rate = 0.004 *24 #/day

# Calculation of the effective plasma concentration of semaglutide for each injection dose
injection_0_25_mg = 0.25 * bioavailability / distribution_volume; #mg/dL
injection_0_5_mg = 0.5 * bioavailability / distribution_volume; #mg/dL
injection_1_mg = 1.0 * bioavailability / distribution_volume; #mg/dL
injection_2_mg = 2.0 * bioavailability / distribution_volume; #mg/dL
injection_1_7_mg = 1.7 * bioavailability / distribution_volume; #mg/dL
injection_2_4_mg = 2.4 * bioavailability / distribution_volume; #mg/dL

# Define the semaglutide doses (mg/dL) (Here, wegovy including dose escalation phase)
injection_doses =[injection_0_25_mg, 
                injection_0_5_mg, 
                injection_1_mg,
                injection_1_7_mg, 
                injection_2_4_mg];

# Define the semaglutide injection times for each dose (days) 
# Each inner array lists the days on which the corresponsing dose is injected
release_times = [[7.0*i for i in 1:4],
                [7.0*i for i in 5:8],
                [7.0*i for i in 9:12],
                [7.0*i for i in 13:16],
                [7.0*i for i in 17:68]];


""" Define reference values, inital values and adjustment times for all variables.
To simulate steady-state, the initial and reference values are set equal"""

# Parameters for body weight
# assume body length of 1.80m
ref_body_weight = initial_body_weight = 107.2 #kg in study 1 from Table 1. BMI > 30
body_weight_adj_time            = 300.0; #days 

# Parameters for food_ingestion
ref_food_ingestion = initial_food_ingestion = 3327.0; #kcal
food_ingestion_adj_time         = 100.0; #days

# Parameters for blood glucose
ref_blood_glucose = initial_blood_glucose  = 110.0; #mg/dL 
blood_glucose_adj_time          = 10.0; #days

# Parameters for GLP-1
ref_glp_1 = initial_glp_1       = 1.088e-6; #mg/dl 
glp_1_adj_time                  = 10.0; #days

# Parameters for semaglutide
ref_semaglutide_subcutaneous = initial_semaglutide_subcutaneous= 0.0; #mg/dl
ref_semaglutide_plasma = initial_semaglutide_plasma = 0.0; #mg/dl

# Parameters for insulin
ref_insulin = initial_insulin   = 7.07e-6 *1.5; #mg/dL  
insulin_adj_time                = 20.0; #days

# Parameters for insulin sensivitiy 
ref_insulin_sen = initial_insulin_sen= 0.5; #dimensionless 
insulin_sen_adj_time            = 50.0; #days

# Parameters for Leptin
ref_leptin = initial_leptin     = 4.0e-3; #mg/dl 
leptin_adj_time                 = 10.0; #days

# Parameters for leptin sensitivity
ref_leptin_sen = initial_leptin_sen = 0.5 #dimensionless 
leptin_sen_adj_time     = 50.0; #days

# Parameters for glucotoxicity
ref_glucotoxicity = initial_glucotoxicity = 0.0; #dimensionless 
glucotoxicity_adj_time          = 30.0; #days

# Parameters for lipotoxicity
ref_lipotoxicity = initial_lipotoxicity = 0.4; #dimensionless 
lipotoxicity_adj_time           = 30.0; #days

# Parameters for beta cell functioning
ref_beta_cell_functioning = initial_beta_cell_functioning = 0.8 #dimensionless
beta_cell_functioning_adj_time          = 50.0; #days

# Parameters for AgRP neuron activity
ref_agrp = initial_agrp         = 1.0; #dimensionless
agrp_adj_time                   = 50.0; #days

# Parameters for POMC neuron activity
ref_pomc = initial_pomc         = 1.0; #dimensionless
pomc_adj_time                   = 15.0; #days

# Parameters for dopamine neuron activity
ref_dopamine = initial_dopamine     = 0.8; #dimensionless
dopamine_adj_time                   = 5.0; #days