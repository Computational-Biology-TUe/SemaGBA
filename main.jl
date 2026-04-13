"""
    <SemaGBA: A system dynamics model of the Semaglutide-responsive Gut-Brain Axis. 
    A model of how the brain and semaglutide regulate appetite and weight.>
    Copyright (C) 2026  Vivan Kennis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

using Pkg, CairoMakie
include("src/Functions.jl") 

Pkg.activate(".")
Pkg.instantiate()

function main(simulations)
    if isempty(simulations)
        simulations = [2,3,4,5,6]
    else
        simulations = parse.(Int, simulations)
    end

    if 2 in simulations
        # run simulation 2
        # Define input and ODE file
        inputfile_diabetes = "../input/Inputs_diabetes.jl"
        inputfile_obesity = "../input/Inputs_obesity.jl" 
        odefile = "../src/ODE_reduced.jl" 

        # Compute normalized changes for diabetic and obesity scenario
        vars_diabetic = normalized_change(inputfile_diabetes, odefile, callback_type = :semaglutide);
        vars_obesity = normalized_change(inputfile_obesity, odefile, callback_type = :semaglutide);

        # Barplot setup
        var_names = ["Body weight*", "Food ingestion*", "Blood glucose*", 
                    "GLP-1*", "Insulin*", "Leptin*", 
                    "Insulin sensitivity", "Leptin sensitivity", "β-cell function", 
                    "Glucotoxicity", "Lipotoxicity"];

        n = length(var_names); # number of variables 
        k=2;                 # number of groups
        cm = 96 / 2.54;      # convert cm to pixels (for figure sizing)

        ys_flat = vcat(vars_diabetic,vars_obesity);  # flatten values for plotting
        xs_flat = repeat(1:n, outer=k);              # x-positions for each varible
        groups = repeat(1:k,inner=n);              # group index for positioning

        # Define colors and patterns for the bars
        colors = [Makie.wong_colors()[1], Makie.wong_colors()[2]];
        patterns=let 
            patternsymbols = ["|", "/"];
            [Pattern(pat, background_color=col, linecolor=:black) for (col,pat) in zip(colors, patternsymbols)];
        end
        elements = [PolyElement(pattern=p) for p in patterns];

        # Create figure

        figure = let f = Figure(size=(17cm,12cm))

            a = Axis(fig[1,1], 
                xticks=(1:n, var_names),        
                xticklabelrotation=pi/6,        
                ylabel="Normalized change from baseline [-]", 
                title="Effect of semaglutide on variables",
                titlesize=18, xlabelsize=16, ylabelsize=16)

            uniquegroups = unique(groups)
            
            dodge_gap = 0.03 
            gap = 0.2 
            width = 1.0 	
            widthplot = width / (length(uniquegroups) + 0.2) 
            
            xs_explicit_flat = compute_x(xs_flat, groups; width, gap, dodge_gap)

            bars = []
            for (g,pat) in zip(uniquegroups, patterns)
                indices = findall(==(g), groups)
                ys_g = ys_flat[indices]
                xs_g = xs_explicit_flat[indices]
                push!(bars, barplot!(a, xs_g, ys_g; color=pat, dodge_gap=dodge_gap, gap=gap, width=widthplot))
            end

            Legend(fig[1,2], bars, ["Diabetes", "Obesity"]) 

            Label(fig[2,1:2], "*Normalized with reference value", fontsize = 12, halign=:right)
            
            # save figure
            figdir = joinpath(@__DIR__, "..", "figures")
            mkpath(figdir)
            save(joinpath(figdir, "fig2_Barplot.png"), f)

            f
        end

    end

    if 3 in simulations
        # run simulation 3
        # Define input and ODE file
        inputfile = "input/Inputs_healthy.jl"
        odefile = "src/ODE_reduced.jl" 

        # Run overeating simulation with (sol_treated) and without (sol_untreated) treatment
        sol_untreated = run_simulation(inputfile, odefile, callback_type= :overeating); 
        sol_treated = run_simulation(inputfile, odefile, callback_type= :overeating_with_semaglutide); 

        extract_values3(sol) = (
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

        variables_untreated = extract_values3(sol_untreated)
        variables_treated = extract_values3(sol_treated)

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
            
            # save figure
            figdir = joinpath(@__DIR__, "..", "figures")
            mkpath(figdir)
            save(joinpath(figdir, "fig3_Overeating.png"), f)   
            f
        end

    end

    if 4 in simulations
        # run simulation 4
        # Define input and ODE file
        inputfile_diabetes = "input/Inputs_diabetes.jl"
        inputfile_obesity = "input/Inputs_obesity.jl" 
        odefile = "src/ODE_extended.jl" 

        # Run diabetes and obesity simulation
        sol_diabetes = run_simulation(inputfile_diabetes, odefile, callback_type= :semaglutide); 
        sol_obesity = run_simulation(inputfile_obesity, odefile, callback_type= :semaglutide);

        # Extract variables
        extract_neurons(sol) = (
            agrp      = [round(u[14], digits=3) for u in sol.u],
            pomc      = [round(u[15], digits=3) for u in sol.u],
            dopamine  = [round(u[16], digits=3) for u in sol.u]);

        neurons_diabetic = extract_neurons(sol_diabetes);
        neurons_obesity  = extract_neurons(sol_obesity);

        cm = 96 / 2.54;      # convert cm to pixels (for figure sizing)

        # Define data, labels, colors and linestyles for the plot
        y_data_diabetic = [neurons_diabetic.agrp, neurons_diabetic.pomc, neurons_diabetic.dopamine];
        y_data_obese = [neurons_obesity.agrp, neurons_obesity.pomc, neurons_obesity.dopamine];
        labels = ["AgRP", "POMC", "Dopamine"];
        colors = [Makie.wong_colors()[1], Makie.wong_colors()[2], Makie.wong_colors()[3]];
        linestyles = [(:dash, :dense), :solid, (:dot, :dense)];
            
        figure = let f = Figure(size=(17cm,14cm), fontsize=12)
            axs = Axis[]
                
            ax1 = Axis(f[1,1], xlabel="Time [days]", ylabel="Neuron activity [-]", title="Diabetic", titlesize =18, xlabelsize = 16, ylabelsize = 16)
            for i in 1:3
                lines!(ax1, sol_diabetes.t, y_data_diabetic[i], color=colors[i], label=labels[i], linestyle=linestyles[i], linewidth = 3) 

            end

            ax2 = Axis(f[2,1], xlabel="Time [days]", ylabel="Neuron activity[-]", title="Obesity", titlesize =18, xlabelsize = 16, ylabelsize = 16)
            for i in 1:3
                lines!(ax2, sol_obesity.t, y_data_obese[i], color=colors[i], label=labels[i], linestyle=linestyles[i], linewidth = 3) 
            end

            Legend(f[1, 2], ax1, orientation=:vertical)
            
            # save figure
            figdir = joinpath(@__DIR__, "..", "figures")
            mkpath(figdir)
            save(joinpath(figdir, "fig4b_Neuron_activity.png"), f)
    
            f     
        end

    end

    if 5 in simulations
        # run simulation 5
        # Define input and ODE file
        inputfile = "input/Inputs_diabetes.jl"
        odefile_reduced = "src/ODE_reduced.jl" 
        odefile_extended = "src/ODE_extended.jl" 

        # Run overeating simulation with (sol_treated) and without (sol_untreated) treatment
        sol_reduced = run_simulation(inputfile, odefile_reduced, callback_type= :semaglutide); 
        sol_extended = run_simulation(inputfile, odefile_extended, callback_type= :semaglutide); 

        extract_values5(sol) = (
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

        variables_reduced = extract_values5(sol_reduced)
        variables_extended = extract_values5(sol_extended)

        cm = 96 / 2.54;      # convert cm to pixels (for figure sizing)

        # Define data, labels, and linestyles for the plot
        y_data_reduced = [variables_reduced.body_weight, variables_reduced.food_ingestion, 
                variables_reduced.blood_glucose, variables_reduced.insulin, variables_reduced.insulin_sen, 
                variables_reduced.glucotoxicity, variables_reduced.leptin, variables_reduced.leptin_sen, 
                variables_reduced.lipotoxicity, variables_reduced.glp1, variables_reduced.semaglu_sub, 
                variables_reduced.beta]

        y_data_extended = [variables_extended.body_weight, variables_extended.food_ingestion, 
                variables_extended.blood_glucose, variables_extended.insulin, variables_extended.insulin_sen, 
                variables_extended.glucotoxicity, variables_extended.leptin, variables_extended.leptin_sen, 
                variables_extended.lipotoxicity, variables_extended.glp1, variables_extended.semaglu_sub, 
                variables_extended.beta]

        labels = ["Body weight\n[kg]", "Food ingestion\n[kcal]", 
                    "Blood glucose\n[mg/dL]", "Insulin\n[mg/dL]", "Insulin sensitivity\n[-]", 
                    "Glucotoxicity\n[-]", "Leptin\n[mg/dL]", "Leptin sensitivity\n[mg/dL]", 
                    "Lipotoxicity\n[-]", "GLP-1\n[mg/dL]", "Semaglutide\n[mg/dL]", 
                    "β-cell functioning\n[mg/dL]"]

        # Create figure
        figure = let f = Figure(size=(15cm,17cm), fontsize=12)
            axs = Axis[]
            tr_text = split("ABCDEFGHIJKL", "")
                for i in eachindex(y_data_reduced)

                    x = (i-1) ÷ 3 
                    y = (i-1) % 3
                    ax = Axis(f[x,y], xlabel="Time [days]", ylabel=labels[i])
                    lines!(ax, sol_reduced.t, y_data_reduced[i], color=Makie.wong_colors()[2], label="Original")
                    lines!(ax, sol_extended.t, y_data_extended[i], color=Makie.wong_colors()[1], label="With extended\n activity", linestyle= :dash)
                    push!(axs, ax)
                    Label(f[x,y, TopLeft()], tr_text[i],padding = (0, 15, 15, 0), font=:bold, fontsize=16)
                end
                Legend(f[4,0:2], axs[1], orientation=:horizontal)
                Label(f[begin-1, 0:2],"Diabetes treatement with 0.5mg semaglutide", fontsize=20, font = :bold)
                
                # save figure
                figdir = joinpath(@__DIR__, "..", "figures")
                mkpath(figdir)
                save(joinpath(figdir, "fig5_Trajectories_diabetic.png"), f)

                f
        end

    end

    if 6 in simulations
        # run simulation 6
        # Define input and ODE file

        inputfile = joinpath(@__DIR__, "input", "Inputs_obesity.jl")
        odefile_reduced = joinpath(@__DIR__, "src", "ODE_reduced.jl")
        odefile_extended = joinpath(@__DIR__, "src", "ODE_extended.jl")

        # Run overeating simulation with (sol_treated) and without (sol_untreated) treatment
        sol_reduced = run_simulation(inputfile, odefile_reduced, callback_type= :semaglutide); 
        sol_extended = run_simulation(inputfile, odefile_extended, callback_type= :semaglutide); 

        extract_values6(sol) = (
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

        variables_reduced = extract_values6(sol_reduced)
        variables_extended = extract_values6(sol_extended)

        cm = 96 / 2.54;      # convert cm to pixels (for figure sizing)

        # Define data, labels, and linestyles for the plot
        y_data_reduced = [variables_reduced.body_weight, variables_reduced.food_ingestion, 
                variables_reduced.blood_glucose, variables_reduced.insulin, variables_reduced.insulin_sen, 
                variables_reduced.glucotoxicity, variables_reduced.leptin, variables_reduced.leptin_sen, 
                variables_reduced.lipotoxicity, variables_reduced.glp1, variables_reduced.semaglu_sub, 
                variables_reduced.beta]

        y_data_extended = [variables_extended.body_weight, variables_extended.food_ingestion, 
                variables_extended.blood_glucose, variables_extended.insulin, variables_extended.insulin_sen, 
                variables_extended.glucotoxicity, variables_extended.leptin, variables_extended.leptin_sen, 
                variables_extended.lipotoxicity, variables_extended.glp1, variables_extended.semaglu_sub, 
                variables_extended.beta]

        labels = ["Body weight\n[kg]", "Food ingestion\n[kcal]", 
                    "Blood glucose\n[mg/dL]", "Insulin\n[mg/dL]", "Insulin sensitivity\n[-]", 
                    "Glucotoxicity\n[-]", "Leptin\n[mg/dL]", "Leptin sensitivity\n[mg/dL]", 
                    "Lipotoxicity\n[-]", "GLP-1\n[mg/dL]", "Semaglutide\n[mg/dL]", 
                    "β-cell functioning\n[mg/dL]"]

        # Create figure
        figure = let f = Figure(size=(15cm,17cm), fontsize=12)
            axs = Axis[]
            tr_text = split("ABCDEFGHIJKL", "")
                for i in eachindex(y_data_reduced)

                    x = (i-1) ÷ 3 
                    y = (i-1) % 3
                    ax = Axis(f[x,y], xlabel="Time [days]", ylabel=labels[i])
                    lines!(ax, sol_reduced.t, y_data_reduced[i], color=Makie.wong_colors()[2], label="Original")
                    lines!(ax, sol_extended.t, y_data_extended[i], color=Makie.wong_colors()[1], label="With extended\n activity", linestyle= :dash)
                    push!(axs, ax)
                    Label(f[x,y, TopLeft()], tr_text[i],padding = (0, 15, 15, 0), font=:bold, fontsize=16)
                end
                Legend(f[4,0:2], axs[1], orientation=:horizontal)
                Label(f[begin-1, 0:2],"Obesity treatement with 2.4mg semaglutide", fontsize=20, font = :bold)
                
                # save figure
                figdir = joinpath(@__DIR__, "..", "figures")
                mkpath(figdir)
                save(joinpath(figdir, "fig66_Trajectories_obesity.png"), f)

                f
        end

    end
    

end


main(["6"])

