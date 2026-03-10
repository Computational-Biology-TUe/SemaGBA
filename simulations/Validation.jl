using CairoMakie

include("../src/Functions.jl") 

# Define input and ODE file
inputfile = "../input/Inputs_healthy.jl"

odefile_reduced = "../src/ODE_reduced.jl" 

include(inputfile) 


# Run simulation with the reduced and extended ODE model
#sol_reduced = run_simulation(inputfile, odefile_reduced, callback_type= :semaglutide); 

extract_values(sol) = (
    body_weight      = [round(u[1], digits=2) for u in sol.u],
    net_energy_intake   = [round(u[2], digits=2) for u in sol.u],
    blood_glucose    = [round(u[3], digits=2) for u in sol.u],
    glp1             = [round(u[4], digits=9) for u in sol.u],
    insulin          = [round(u[5], digits=7) for u in sol.u],
    insulin_sen      = [round(u[6], digits=3) for u in sol.u],
    leptin           = [round(u[7], digits=6) for u in sol.u],
    leptin_sen       = [round(u[8], digits=3) for u in sol.u],
    glucotoxicity    = [round(u[9], digits=3) for u in sol.u],
    lipotoxicity     = [round(u[10], digits=3) for u in sol.u],
    semaglu_sub      = [round(u[11], digits=9) for u in sol.u],
    semaglu_plasma   = [round(u[12], digits=9) for u in sol.u],
    beta             = [round(u[13], digits=3) for u in sol.u])

#variables_reduced = extract_values(sol_reduced)

cm = 96 / 2.54;      # convert cm to pixels (for figure sizing)

# Define data, labels, and linestyles for the plot
#= y_data_reduced = [variables_reduced.body_weight, variables_reduced.net_energy_intake, 
        variables_reduced.blood_glucose, variables_reduced.insulin, variables_reduced.insulin_sen, 
        variables_reduced.glucotoxicity, variables_reduced.leptin, variables_reduced.leptin_sen, 
        variables_reduced.lipotoxicity, variables_reduced.glp1, variables_reduced.semaglu_plasma, 
        variables_reduced.beta]

labels = ["Body weight\n[kg]", "Net energy intake\n[kcal]", 
            "Blood glucose\n[mg/dL]", "Insulin\n[μU/mL]", "Insulin sensitivity\n[-]", 
            "Glucotoxicity\n[-]", "Leptin\n[mg/dL]", "Leptin sensitivity\n[-]", 
            "Lipotoxicity\n[-]", "GLP-1\n[mg/dL]", "Semaglutide\n[mg/dL]", 
            "β-cell functioning\n[-]"]
 =#
# Create figure
#= figure = let f = Figure(size=(15cm,17cm), fontsize=12)
    axs = Axis[]
    tr_text = split("ABCDEFGHIJKL", "")
        for i in eachindex(y_data_reduced)

            x = (i-1) ÷ 3 
            y = (i-1) % 3
            ax = Axis(f[x,y], xlabel="Time [days]", ylabel=labels[i])
            lines!(ax, sol_reduced.t, y_data_reduced[i], color=Makie.wong_colors()[2], label="Reduced model")
            push!(axs, ax)
            Label(f[x,y, TopLeft()], tr_text[i],padding = (0, 15, 15, 0), font=:bold, fontsize=16)
        end
        Legend(f[4,0:2], axs[1], orientation=:horizontal)
        Label(f[begin-1, 0:2],"Validation of model with extremely high body weight", fontsize=18, font = :bold)
        
        # save figure
        figdir = joinpath(@__DIR__, "..", "figures")
        mkpath(figdir)
        save(joinpath(figdir, "fig9_validation_extreme.png"), f)

        f
end =#

sol_1 = run_simulation(inputfile, odefile_reduced, callback_type= :semaglutide, 
        release_times = [[7]], 
        injection_doses = [injection_0_25_mg])

sol_05 = run_simulation(inputfile, odefile_reduced, callback_type= :semaglutide,
        release_times = [[7.0*i for i in 1:4],[7.0*i for i in 5:31]], 
        injection_doses =[injection_0_25_mg, injection_0_5_mg])

sol_10 = run_simulation(inputfile, odefile_reduced, callback_type= :semaglutide,
        release_times = [[7.0*i for i in 1:4],[7.0*i for i in 5:8],[7 *i for i in 9:31]], 
        injection_doses =[injection_0_25_mg, injection_0_5_mg, injection_1_mg])

sol_oz = run_simulation(inputfile, odefile_reduced, callback_type= :semaglutide,
        release_times = [[7.0*i for i in 1:4],[7.0*i for i in 5:8],[7 *i for i in 9:31]], 
        injection_doses =[injection_0_25_mg, injection_0_5_mg, injection_1_mg])

sol_weg = run_simulation(inputfile, odefile_reduced, callback_type= :semaglutide,
        release_times = [[7.0*i for i in 1:4],[7.0*i for i in 5:8],[7.0*i for i in 9:12],[7.0*i for i in 13:16],[7.0*i for i in 17:31]], 
        injection_doses =[injection_0_25_mg, injection_0_5_mg, injection_1_mg, injection_1_7_mg, injection_2_4_mg])

variables_1 = extract_values(sol_1)
variables_05 = extract_values(sol_05)
variables_10 = extract_values(sol_10)

variables_oz = extract_values(sol_oz)
variables_weg = extract_values(sol_weg)


# Create figure 
figure = let f = Figure(fontsize=12)    
    ax = Axis(f[1, 1], xlabel = "Time [days]", ylabel = "Semaglutide [mg/dL]", title="Semaglutide concentration in plasma", titlesize =18,  xlabelsize = 16, ylabelsize = 16)
    lines!(ax, sol_weg.t, variables_weg.semaglu_plasma, color = Makie.wong_colors()[3], label = "Wegovy dosing", linewidth = 3)
    lines!(ax, sol_oz.t, variables_oz.semaglu_plasma, color = Makie.wong_colors()[1], label = "Ozempic dosing", linewidth = 2)
    lines!(ax, sol_1.t, variables_1.semaglu_plasma, color = Makie.wong_colors()[2], label = "One injection", linewidth = 2)


    Legend(f[2, 1], ax, orientation = :horizontal)

    # Save figure
    figdir = joinpath(@__DIR__, "..", "figures")
    mkpath(figdir)
    save(joinpath(figdir, "fig10_semaglutide_plasma.png"), f)
    f
end

# maximale concentratie
max_plasma_conc = maximum(variables_1.semaglu_plasma)
idx = argmax(variables_1.semaglu_plasma)
time_max = sol_1.t[idx]
println("Max plasma concentration = ", max_plasma_conc)
println("Reached at time = ", time_max, " days")

# Percentage based on simulated plasma values (one 0.5 mg injection at day 7)
injection_time = 7.0
target_time = injection_time + 7 * 7
idx_target = argmin(abs.(sol_1.t .- target_time))

# Use unrounded state values for the concentration comparison
plasma_series = [u[12] for u in sol_1.u]
plasma_at_target = plasma_series[idx_target]
percent_remaining_vs_injected = 100 * plasma_at_target / injection_0_5_mg

println("Referentie (injection_0_5_mg): ", injection_0_5_mg, " mg/dL")
println("Plasma op dag ", sol_1.t[idx_target], " = ", plasma_at_target, " mg/dL")
println("Percentage over na 7 weken na injectie: ", percent_remaining_vs_injected, " %")

# Filter time points for 70 < t < 30*7 days
t_min = 7*16
t_max = 28 * 7.0

figure = let f = Figure(fontsize=12)    
    ax = Axis(f[1, 1], xlabel = "Time [days]", ylabel = "Semaglutide [mg/dL]", 
              title="Semaglutide concentration in plasma", titlesize=18, 
              xlabelsize=16, ylabelsize=16)
    
        # Filter indices in the requested time window
        idx_weg = (sol_10.t .> t_min) .& (sol_10.t .< t_max)
    
    lines!(ax, sol_10.t[idx_weg], variables_10.semaglu_plasma[idx_weg], 
           color=Makie.wong_colors()[3], label="Wegovy dosing", linewidth=3)
    
    Legend(f[2, 1], ax, orientation=:horizontal)

    f
end
