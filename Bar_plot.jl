using CairoMakie
include("Functions.jl") 

# Define input and ODE file
inputfile_diabetes = "Inputs_diabetes.jl"
inputfile_obesity = "Inputs_obesity.jl" 
odefile = "ODE_simplified.jl" 

# Compute normalized changes for diabetic and obesity scenario
vars_diabetic = normalized_change(inputfile_diabetes, odefile, callback_type = "semaglutide");
vars_obesity = normalized_change(inputfile_obesity, odefile, callback_type = "semaglutide");

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
let 
    fig = Figure(size=(17cm,12cm))

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
    save("barplot.png", fig)

	fig
end
