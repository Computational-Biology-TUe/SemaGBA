using CairoMakie
include("../src/Functions.jl") 

# Define input and ODE file
inputfile_diabetes = "../input/Inputs_diabetes.jl"
inputfile_obesity = "../input/Inputs_obesity.jl" 
odefile = "../src/ODE_neural.jl" 

function compute_effects(variables, inputfile)
        include(inputfile)

        # Effects on AgRP neuron activity 
        list_effect_i_on_agrp = [
                do_interpolate(effect_i_on_agrp_interp, (point / ref_insulin), xs_i_agrp)
                for point in variables.insulin
            ]
        list_effect_l_on_agrp = [
                do_interpolate(effect_l_on_agrp_interp, (point / ref_leptin), xs_l_agrp)
                for point in variables.leptin
            ]
        list_effect_do_on_agrp = [
                do_interpolate(effect_do_on_agrp_interp, (point - ref_dopamine), xs_do_agrp)
                for point in variables.dopamine
            ]

        # Effects on POMC neuron activity 
        list_effect_i_on_pomc = [
                do_interpolate(effect_i_on_pomc_interp, (point / ref_insulin), xs_i_pomc)
                for point in variables.insulin
            ]
        list_effect_l_on_pomc = [
                do_interpolate(effect_l_on_pomc_interp, (point / ref_leptin), xs_l_pomc)
                for point in variables.leptin
            ]
        list_effect_glp_on_pomc = [
                do_interpolate(effect_glp_on_pomc_interp, ((point+sema) / ref_glp_1), xs_glp_pomc)
                for (point, sema) in zip(variables.glp1, variables.semaglu_plasma)
            ]

        list_effect_do_on_pomc = [
                do_interpolate(effect_do_on_pomc_interp, (point - ref_dopamine), xs_do_pomc)
                for point in variables.dopamine
            ]

        list_effect_agrp_on_pomc = [
                do_interpolate(effect_agrp_on_pomc_interp, (point - ref_agrp), xs_agrp_pomc)
                for point in variables.agrp
            ]

        # Effects on dopamine neuron activity 
        list_effect_i_on_dopamine = [
                do_interpolate(effect_i_on_do_interp, (point / ref_insulin), xs_i_do)
                for point in variables.insulin
            ]
        list_effect_l_on_dopamine = [
                do_interpolate(effect_l_on_do_interp, (point / ref_leptin), xs_l_do)
                for point in variables.leptin
            ]
        list_effect_glp_on_dopamine = [
                do_interpolate(effect_glp_on_do_interp, ((point+sema) / ref_glp_1), xs_glp_do)
                for (point, sema) in zip(variables.glp1, variables.semaglu_plasma)
            ]
         # Define effects and labels for plotting
        effects_agrp = [list_effect_i_on_agrp, list_effect_l_on_agrp, list_effect_do_on_agrp]
        labels_agrp = ["insulin", "leptin", "dopamine\nneuron activity"]

        effects_pomc= [list_effect_i_on_pomc, list_effect_l_on_pomc, list_effect_glp_on_pomc, list_effect_do_on_pomc, list_effect_agrp_on_pomc]
        labels_pomc = ["insulin", "leptin", "GLP-1", "dopamine\nneuron activity", "AgRP\nneuron activity"]

        effects_dopamine = [list_effect_i_on_dopamine, list_effect_l_on_dopamine, list_effect_glp_on_dopamine]
        labels_dopamine = ["insulin", "leptin", "GLP-1"]

        effects = [effects_agrp, effects_pomc, effects_dopamine]
        labels = [labels_agrp, labels_pomc, labels_dopamine]

        neurons = [variables.agrp, variables.pomc, variables.dopamine]
        return effects, labels, neurons
    end

function effect_plot(inputfile_diabetes, inputfile_obesity, odefile)
    # Run simulation
    sol_diabetes = run_simulation(inputfile_diabetes, odefile, callback_type= :semaglutide);
    sol_obesity = run_simulation(inputfile_obesity, odefile, callback_type= :semaglutide);

    # Extract values of all variables for diabetes and obesity simulation
    extract_values(sol) = (
        glp1             = [u[4] for u in sol.u],
        insulin          = [u[5] for u in sol.u],
        leptin           = [u[7] for u in sol.u],
        semaglu_plasma   = [u[12] for u in sol.u],
        agrp             = [u[14] for u in sol.u],
        pomc             = [u[15] for u in sol.u],
        dopamine         = [u[16] for u in sol.u]);

    variables_diabetes = extract_values(sol_diabetes)
    variables_obesity = extract_values(sol_obesity)
    
    effects_diabetes, labels, neurons_diabetes = compute_effects(variables_diabetes, inputfile_diabetes)
    effects_obesity, _, neurons_obesity = compute_effects(variables_obesity, inputfile_obesity)
    
    # Define labels and color map for plotting

    color_map = Dict(
        "insulin" => Makie.wong_colors()[1],
        "leptin" => Makie.wong_colors()[2],
        "GLP-1" => Makie.wong_colors()[3],
        "dopamine\nneuron activity" => Makie.wong_colors()[4],
        "AgRP\nneuron activity" => Makie.wong_colors()[7])
    
    # Combine data for diabetes and obesity simulation for plotting
    y2_labels = ["Resulting AgRP neuron activity", "Resulting POMC neuron activity", "Resulting dopamine neuron activity"]
    sols = [sol_diabetes, sol_obesity]
    effects_all = [effects_diabetes, effects_obesity]
    neurons_all = [neurons_diabetes, neurons_obesity]
    title = ["Diabetes", "Obesity"]
    
    # Create figure 
    figure = let f = Figure(fontsize=10, size=(15cm,17cm))
        all_plots = Dict{String, Any}() 
        tr_text = split("ABCDEF", "")

        Label(f[0, 1:2], "Contributions to changes in neuronal activity over time", fontsize = 18, font = :bold)
        
        # Create three subplots for AgRP, POMC and dopamine neuron activity
        for col in 1:2
            for row in 1:3
                ax = Axis(f[row, col], xlabel = "Time [days]", ylabel = "Effect [-]")
                for effect in eachindex(effects_all[col][row])
                    p = lines!(ax, sols[col].t, effects_all[col][row][effect], color = color_map[labels[row][effect]], label = labels[row][effect])
                    if !haskey(all_plots, labels[row][effect])
                        all_plots[labels[row][effect]] = p
                    end
                end

                # Plot of resulting neuronal activity
                if row <= 1
                    ax2 = Axis(f[row, col], yaxisposition = :right, ylabel = y2_labels[row], title=title[col], titlesize=14)
                else
                    ax2 = Axis(f[row, col], yaxisposition = :right, ylabel = y2_labels[row])
                end
                hidespines!(ax2)
                hidexdecorations!(ax2)
                p2 = lines!(ax2, sols[col].t, neurons_all[col][row], color = :black, linestyle = :dash)
                if !haskey(all_plots, "Resulting neuron activity")
                    all_plots["Resulting neuron activity"] = p2
                end
                Label(f[row, col, TopLeft()], tr_text[row + (col-1)*3], padding = (0, 15, 15, 0), font=:bold, fontsize=14)

            end         
        end
            # Add legend
            legend_labels = filter(x -> x != "Resulting neuron activity", collect(keys(all_plots)))
            legend_labels = vcat(legend_labels, "Resulting neuron activity")
            legend_plots = [all_plots[label] for label in legend_labels]
            Legend(f[4, 1:2], legend_plots, legend_labels, orientation = :horizontal, labelsize=9, nbanks = 1)
            
            # Save figure
            figdir = joinpath(@__DIR__, "..", "figures")
            mkpath(figdir)
            save(joinpath(figdir, "fig8_contribution_neuron.png"), f)
        f
    end
end

effect_plot(inputfile_diabetes, inputfile_obesity, odefile)
