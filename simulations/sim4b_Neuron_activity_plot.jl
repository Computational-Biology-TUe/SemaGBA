using CairoMakie
include("../src/Functions.jl") 

# Define input and ODE file
inputfile_diabetes = "../input/Inputs_diabetes.jl"
inputfile_obesity = "../input/Inputs_obesity.jl" 
odefile = "../src/ODE_extended.jl" 

# # Run simulation with the diabetes and obesity inputs
sol_diabetes = run_simulation(inputfile_diabetes, odefile, callback_type= :semaglutide); 
sol_obesity = run_simulation(inputfile_obesity, odefile, callback_type= :semaglutide);

# Extract variables
extract_neurons(sol) = (
    agrp      = [round(u[14], digits=3) for u in sol.u],
    pomc      = [round(u[15], digits=3) for u in sol.u],
    dopamine  = [round(u[16], digits=3) for u in sol.u]);

neurons_diabetes = extract_neurons(sol_diabetes);
neurons_obesity  = extract_neurons(sol_obesity);

cm = 96 / 2.54;      # convert cm to pixels (for figure sizing)

# Define data, labels, colors and linestyles for the plot
y_data_diabetes = [neurons_diabetes.agrp, neurons_diabetes.pomc, neurons_diabetes.dopamine];
y_data_obesity = [neurons_obesity.agrp, neurons_obesity.pomc, neurons_obesity.dopamine];
labels = ["AgRP", "POMC", "Dopamine"];
colors = [Makie.wong_colors()[1], Makie.wong_colors()[2], Makie.wong_colors()[3]];
linestyles = [(:dash, :dense), :solid, (:dot, :dense)];
     
figure = let f = Figure(size=(17cm,14cm), fontsize=12)
    axs = Axis[]
    
    # Plot diabetes neuron activity
    ax1 = Axis(f[1,1], xlabel="Time [days]", ylabel="Neuron activity [-]", title="Diabetes", titlesize =18, xlabelsize = 16, ylabelsize = 16)
    for i in 1:3
        lines!(ax1, sol_diabetes.t, y_data_diabetes[i], color=colors[i], label=labels[i], linestyle=linestyles[i], linewidth = 3) 

    end
    
    # Plot obesity neuron activity
    ax2 = Axis(f[2,1], xlabel="Time [days]", ylabel="Neuron activity[-]", title="Obesity", titlesize =18, xlabelsize = 16, ylabelsize = 16)
    for i in 1:3
        lines!(ax2, sol_obesity.t, y_data_obesity[i], color=colors[i], label=labels[i], linestyle=linestyles[i], linewidth = 3) 
    end

    Legend(f[1, 2], ax1, orientation=:vertical)

    # save figure
    figdir = joinpath(@__DIR__, "..", "figures")
    mkpath(figdir)
    save(joinpath(figdir, "fig4b_Neuron_activity.png"), f)
    
    f
        
end
