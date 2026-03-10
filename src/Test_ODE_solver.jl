inputfile_diabetes = "../input/Inputs_diabetes.jl"
inputfile_obesity = "../input/Inputs_obesity.jl" 
odefile = "../src/ODE_extended.jl" 
callback_type = :semaglutide

include(inputfile_obesity)
include(odefile)
    
# Define ODE
tspan = (0.0, t_days)
u0 = [initial_body_weight, initial_net_energy_intake, initial_blood_glucose, 
        initial_glp_1,  initial_insulin, initial_insulin_sen, 
        initial_leptin, initial_leptin_sen, initial_glucotoxicity, 
        initial_lipotoxicity, initial_semaglutide_subcutaneous, initial_semaglutide_plasma, 
        initial_beta_cell_functioning, initial_agrp, initial_pomc, initial_dopamine];

p = (ref_body_weight, ref_net_energy_intake, ref_blood_glucose,
        ref_glp_1, ref_insulin, ref_insulin_sen, 
        ref_leptin, ref_leptin_sen, ref_glucotoxicity,
        ref_lipotoxicity, ref_semaglutide_subcutaneous, 
        ref_beta_cell_functioning, 
        ref_agrp, ref_pomc,  ref_dopamine,
        absorption_rate, elimination_rate);
    
# Concatenate injection arrays 
injections_time_points = vcat(release_times...)
injections_doses = vcat([fill(dose, length(t)) for (dose,t) in zip(injection_doses, release_times)]...)
    
# Determine callback
if callback_type == :none
    cb = CallbackSet()

elseif callback_type == :semaglutide        
    cb = CallbackSet(semaglutide_callback(injections_time_points, injections_doses))

elseif callback_type == :overeating
    cb = CallbackSet(overeating_callback())

elseif callback_type == :overeating_with_semaglutide
    cb_sema = semaglutide_callback(injections_time_points, injections_doses)
    cb_food = overeating_callback()
    cb = CallbackSet(cb_food,cb_sema);
else 
    error("Unknown type. Use \"none\",\"semaglutide\", \"overeating\"  or \"overeating with semaglutide\".")
end 

# Run simulation
prob = ODEProblem(ODEs!, u0, tspan, p)
sol = solve(prob, alg_hints=:[stiff], tstops=injections_time_points, saveat=0.1, callback=cb);


sol1 = solve(prob, Tsit5(), tstops=injections_time_points, saveat=0.1, callback=cb)
sol2=  solve(prob,Rodas5(), tstops=injections_time_points, saveat=0.1, callback=cb)
