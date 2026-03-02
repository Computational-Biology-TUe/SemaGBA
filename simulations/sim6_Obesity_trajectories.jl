using CairoMakie

include("../src/Functions.jl") 

# Define input and ODE file
inputfile = "../input/Inputs_obesity.jl"
odefile_simplified = "../src/ODE_simplified.jl" 
odefile_neural = "../src/ODE_neural.jl" 

# Run overeating simulation with (sol_treated) and without (sol_untreated) treatment
sol_simplified = run_simulation(inputfile, odefile_simplified, callback_type= :semaglutide); 
sol_neural = run_simulation(inputfile, odefile_neural, callback_type= :semaglutide); 

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

variables_simplified = extract_values(sol_simplified)
variables_neural = extract_values(sol_neural)

cm = 96 / 2.54;      # convert cm to pixels (for figure sizing)

# Define data, labels, and linestyles for the plot
y_data_simplified = [variables_simplified.body_weight, variables_simplified.net_energy_intake, 
        variables_simplified.blood_glucose, variables_simplified.insulin, variables_simplified.insulin_sen, 
        variables_simplified.glucotoxicity, variables_simplified.leptin, variables_simplified.leptin_sen, 
        variables_simplified.lipotoxicity, variables_simplified.glp1, variables_simplified.semaglu_plasma, 
        variables_simplified.beta]

y_data_neural = [variables_neural.body_weight, variables_neural.net_energy_intake, 
        variables_neural.blood_glucose, variables_neural.insulin, variables_neural.insulin_sen, 
        variables_neural.glucotoxicity, variables_neural.leptin, variables_neural.leptin_sen, 
        variables_neural.lipotoxicity, variables_neural.glp1, variables_neural.semaglu_plasma, 
        variables_neural.beta]

labels = ["Body weight\n[kg]", "Net energy intake\n[kcal]", 
            "Blood glucose\n[mg/dL]", "Insulin\n[μU/mL]", "Insulin sensitivity\n[-]", 
            "Glucotoxicity\n[-]", "Leptin\n[mg/dL]", "Leptin sensitivity\n[-]", 
            "Lipotoxicity\n[-]", "GLP-1\n[mg/dL]", "Semaglutide\n[mg/dL]", 
            "β-cell functioning\n[-]"]

# Create figure
figure = let f = Figure(size=(15cm,17cm), fontsize=12)
    axs = Axis[]
    tr_text = split("ABCDEFGHIJKL", "")
        for i in eachindex(y_data_simplified)

            x = (i-1) ÷ 3 
            y = (i-1) % 3
            ax = Axis(f[x,y], xlabel="Time [days]", ylabel=labels[i])
            lines!(ax, sol_simplified.t, y_data_simplified[i], color=Makie.wong_colors()[2], label="Reduced model")
            lines!(ax, sol_neural.t, y_data_neural[i], color=Makie.wong_colors()[1], label="Extended model", linestyle= :dash)
            push!(axs, ax)
            Label(f[x,y, TopLeft()], tr_text[i],padding = (0, 15, 15, 0), font=:bold, fontsize=16)
        end
        Legend(f[4,0:2], axs[1], orientation=:horizontal)
        Label(f[begin-1, 0:2],"Obesity treatement with 2.4mg semaglutide", fontsize=16, font = :bold)
        
        # save figure
        figdir = joinpath(@__DIR__, "..", "figures")
        mkpath(figdir)
        save(joinpath(figdir, "fig6_Trajectories_obesity.png"), f)

        f
end

# Calculate percent difference in net energy intake between simplified and neural model
percent_diff = (variables_simplified.net_energy_intake .- variables_neural.net_energy_intake) ./
                       variables_simplified.net_energy_intake .*100

maximum(abs.(percent_diff))