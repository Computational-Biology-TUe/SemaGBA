using CairoMakie

include("../src/Functions.jl") 

# Define input and ODE file
inputfile = "../input/Inputs_diabetes.jl"

odefile_reduced = "../src/ODE_reduced.jl" 
odefile_extended = "../src/ODE_extended.jl" 

# Run simulation with the reduced and extended ODE model
sol_reduced = run_simulation(inputfile, odefile_reduced, callback_type= :semaglutide); 
sol_extended = run_simulation(inputfile, odefile_extended, callback_type= :semaglutide); 

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

variables_reduced = extract_values(sol_reduced)
variables_extended = extract_values(sol_extended)

cm = 96 / 2.54;      # convert cm to pixels (for figure sizing)

# Define data, labels, and linestyles for the plot
y_data_reduced = [variables_reduced.body_weight, variables_reduced.net_energy_intake, 
        variables_reduced.blood_glucose, variables_reduced.insulin, variables_reduced.insulin_sen, 
        variables_reduced.glucotoxicity, variables_reduced.leptin, variables_reduced.leptin_sen, 
        variables_reduced.lipotoxicity, variables_reduced.glp1, variables_reduced.semaglu_plasma, 
        variables_reduced.beta]

y_data_extended = [variables_extended.body_weight, variables_extended.net_energy_intake, 
        variables_extended.blood_glucose, variables_extended.insulin, variables_extended.insulin_sen, 
        variables_extended.glucotoxicity, variables_extended.leptin, variables_extended.leptin_sen, 
        variables_extended.lipotoxicity, variables_extended.glp1, variables_extended.semaglu_plasma, 
        variables_extended.beta]

labels = ["Body weight\n[kg]", "Net energy intake\n[kcal]", 
            "Blood glucose\n[mg/dL]", "Insulin\n[μU/mL]", "Insulin sensitivity\n[-]", 
            "Glucotoxicity\n[-]", "Leptin\n[mg/dL]", "Leptin sensitivity\n[-]", 
            "Lipotoxicity\n[-]", "GLP-1\n[mg/dL]", "Semaglutide\n[mg/dL]", 
            "β-cell functioning\n[-]"]

# Create figure
figure = let f = Figure(size=(15cm,17cm), fontsize=12)
    axs = Axis[]
    tr_text = split("ABCDEFGHIJKL", "")
        for i in eachindex(y_data_reduced)

            x = (i-1) ÷ 3 
            y = (i-1) % 3
            ax = Axis(f[x,y], xlabel="Time [days]", ylabel=labels[i])
            lines!(ax, sol_reduced.t, y_data_reduced[i], color=Makie.wong_colors()[2], label="Reduced model")
            lines!(ax, sol_extended.t, y_data_extended[i], color=Makie.wong_colors()[1], label="Extended model", linestyle= :dash)
            push!(axs, ax)
            Label(f[x,y, TopLeft()], tr_text[i],padding = (0, 15, 15, 0), font=:bold, fontsize=16)
        end
        Legend(f[4,0:2], axs[1], orientation=:horizontal)
        Label(f[begin-1, 0:2],"Diabetes treatement with 0.5mg semaglutide", fontsize=20, font = :bold)
        
        # save figure
        figdir = joinpath(@__DIR__, "..", "figures")
        mkpath(figdir)
        save(joinpath(figdir, "fig5_Trajectories_diabetes.png"), f)

        f
end

# Calculate percent difference in net energy intake between reduced and extended model
percent_diff = (variables_reduced.net_energy_intake .- variables_extended.net_energy_intake) ./
                       variables_reduced.net_energy_intake .*100

maximum(abs.(percent_diff))


