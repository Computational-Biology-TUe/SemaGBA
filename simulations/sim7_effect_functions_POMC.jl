using Interpolations, Plots, CairoMakie

# Define input and ODE file
include("../src/Data.jl");
include("../src/Functions.jl")

# Define data and labels for the plot
xdata = [xs_i_pomc, xs_l_pomc, xs_glp_pomc, xs_agrp_pomc, xs_do_pomc]
ydata = [ys_i_pomc, ys_l_pomc, ys_glp_pomc, ys_agrp_pomc, ys_do_pomc]
xlabel = ["Current / reference insulin [-]", "Current / reference leptin [-]", "Current / reference GLP-1 [-]", "Current - reference AgRP [-]", "Current - reference dopamine [-]"]
tr_text = split("ABCDEFGHIJKL", "")

# Create figure
figure = let f = Figure(fontsize=12)
    axs = Axis[]

    for i in eachindex(ydata)
        x = (i-1) ÷ 3 
        y = (i-1) % 3
        ax = Axis(f[x,y], xlabel= xlabel[i], ylabel="Effect on POMC activity [-]")
        CairoMakie.lines!(ax, xdata[i], linear_interpolation(xdata[i], ydata[i]).(xdata[i]), color=Makie.wong_colors()[1], label = "Interpolated fit")
        push!(axs, ax)
        Label(f[x,y, TopLeft()], tr_text[i], padding = (0, 15, 15, 0), font=:bold, fontsize=16)
    end
            
    Label(f[begin-1, 0:2],"Effect functions for POMC neuron activity", fontsize=15, font = :bold)

    # save figure
    figdir = joinpath(@__DIR__, "..", "figures")
    mkpath(figdir)
    save(joinpath(figdir, "fig7_Effect_functions_POMC.png"), f)

    f

end
