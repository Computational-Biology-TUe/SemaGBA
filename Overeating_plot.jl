using CairoMakie

# Define input and ODE file
inputfile = "Inputs_healthy.jl"
odefile = "ODE_simplified.jl" 


# Run overeating simulation with (sol_treated) and without (sol_untreated) treatment
sol_untreated = run_simulation(inputfile, odefile, callback_type= "overeating"); 
sol_treated = run_simulation(inputfile, odefile, callback_type= "overeating with semaglutide"); 

extract_values(sol) = (
    body_weight      = [round(u[1], digits=2) for u in sol.u],
    food_ingestion   = [round(u[2], digits=2) for u in sol.u],
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

variables_untreated = extract_values(sol_untreated)
variables_treated = extract_values(sol_treated)

cm = 96 / 2.54;      # convert cm to pixels (for figure sizing)

# Define data, labels, and linestyles for the plot
y_data = [variables_untreated.body_weight, variables_untreated.food_ingestion, 
        variables_untreated.blood_glucose, variables_untreated.insulin, variables_untreated.insulin_sen, 
        variables_untreated.glucotoxicity, variables_untreated.leptin, variables_untreated.leptin_sen, 
        variables_untreated.lipotoxicity, variables_untreated.glp1, variables_untreated.semaglu_sub, 
        variables_untreated.beta]

y_data_treated = [variables_treated.body_weight, variables_treated.food_ingestion, 
        variables_treated.blood_glucose, variables_treated.insulin, variables_treated.insulin_sen, 
        variables_treated.glucotoxicity, variables_treated.leptin, variables_treated.leptin_sen, 
        variables_treated.lipotoxicity, variables_treated.glp1, variables_treated.semaglu_sub, 
        variables_treated.beta]

labels = ["Body weight\n[kg]", "Food ingestion\n[kcal]", 
            "Blood glucose\n[mg/dL]", "Insulin\n[mg/dL]", "Insulin sensitivity\n[-]", 
            "Glucotoxicity\n[-]", "Leptin\n[mg/dL]", "Leptin sensitivity\n[mg/dL]", 
            "Lipotoxicity\n[-]", "GLP-1\n[mg/dL]", "Semaglutide\n[mg/dL]", 
            "β-cell functioning\n[mg/dL]"]

# Create figure
figure = let f = Figure(size=(15cm,17cm), fontsize=12)
    axs = Axis[]

    tr_text = split("ABCDEFGHIJKL", "")
    for i in eachindex(y_data)
        x = (i-1) ÷ 3 
        y = (i-1) % 3
        ax = Axis(f[x,y], xlabel="time [days]", ylabel=labels[i])
        lines!(ax, sol_untreated.t, y_data[i], color=Makie.wong_colors()[2], label="Untreated")
        lines!(ax, sol_treated.t, y_data_treated[i], color=Makie.wong_colors()[1], label="Treated with 1.0mg semaglutide")
        push!(axs, ax)
        Label(f[x,y, TopLeft()], tr_text[i],padding = (0, 15, 15, 0), font=:bold, fontsize=16)
    end
    Label(f[begin-1, 0:2],"Overeating simulation", fontsize=18, font = :bold)
    Legend(f[4,0:2], axs[1], orientation=:horizontal)
    save("overeating.png", f)
    f
end

