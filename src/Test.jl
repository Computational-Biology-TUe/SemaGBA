include("Functions.jl") 
using Plots
sol= sol_simplified
body_weight =   [u[1] for u in sol.u];
food_ingestion = [u[2] for u in sol.u];
blood_glucose =  [u[3] for u in sol.u];
glp1 =      [u[4] for u in sol.u];
insulin =   [u[5] for u in sol.u];
insulin_sen = [u[6] for u in sol.u];
leptin =    [u[7] for u in sol.u];
leptin_sen =    [u[8] for u in sol.u];
glucotoxicity = [u[9] for u in sol.u];
lipotoxicity =  [u[10] for u in sol.u];
semaglu_sub = [u[11] for u in sol.u]
semaglu_plasma = [u[12] for u in sol.u]
beta_cell = [u[13] for u in sol.u]; 

#= list_effect_bg_on_i = [
        do_interpolate(effect_bg_on_i_interp, (point / ref_blood_glucose), xs_bg_i)
        for point in blood_glucose
    ]
list_effect_is_on_i = [
        do_interpolate(effect_is_on_i_interp, (point - ref_insulin_sen), xs_is_i)
        for point in insulin_sen
    ]
    
effect_glp1_over_time = Float64[]

for (i, bg) in enumerate(blood_glucose)
        glp = glp1[i]
        sema = semaglu_plasma[i]
        bg = blood_glucose[i]


        #= if bg == ref_blood_glucose
            effect_val = 1.0
        elseif (bg / 90 > 1.05) && ((glp+sema)/ref_glp_1>1.01)
            effect_val = ((glp + sema) / ref_glp_1) * (bg / 90)
        
        else
            effect_val = 1.0
        end =#

        effect_glp_on_i1 = do_interpolate(effect_glp_on_i_interp,((glp + sema) / ref_glp_1),xs_glp_i);
        glucose_weight = clamp((bg - (90 * 1.00)) / (90 * (1.10 - 1.00)), 0.0, 1.0);
        effect_val = 1.0 + (effect_glp_on_i1-1) * glucose_weight

        push!(effect_glp1_over_time, effect_val)
end


list_effect_b_on_i = [
        do_interpolate(effect_b_on_i_interp, (ref_beta_cell_functioning - point), xs_b_i)

        for point in beta_cell
    ]

Plots.plot(sol.t, list_effect_bg_on_i, label="effect blood glucose on insulin")
Plots.plot!(sol.t, list_effect_is_on_i, label="effect insulin sensitivity on insulin")
Plots.plot!(sol.t, effect_glp1_over_time, label="effect GLP-1 + semaglutide on insulin")
Plots.plot!(sol.t, list_effect_b_on_i, label="effect beta-cell functioning on insulin", legend =:right)

Plots.plot!(twinx(),sol.t, insulin, label="insulin", linestyle=:dash, color=:black, legend=:topleft) =#

list_effect_li_on_is = [
        do_interpolate(effect_li_on_is_interp, (point - ref_lipotoxicity), xs_li_is, ys_li_is)
        for point in lipotoxicity
    ]
list_effect_g_on_is = [
        do_interpolate(effect_g_on_is_interp,(point - ref_glucotoxicity), xs_g_is, ys_g_is)
        for point in glucotoxicity
    ]

Plots.plot(sol.t, list_effect_li_on_is, label="effect lipotoxicity on insulin sensitivity")
Plots.plot!(sol.t, list_effect_g_on_is, label="effect glucotoxicity on insulin sensitivity")
Plots.plot!(twinx(),sol.t, insulin_sen, label="insulin sensitivity", linestyle=:dash, color=:black, legend=:bottom)


using Pkg, DifferentialEquations, Plots, DiffEqCallbacks

variable_name = :dopamine
sol = sol_neural
body_weight =   [u[1] for u in sol.u];
food_ingestion = [u[2] for u in sol.u];
blood_glucose =  [u[3] for u in sol.u];
glp1 =      [u[4] for u in sol.u];
insulin =   [u[5] for u in sol.u];
insulin_sen = [u[6] for u in sol.u];
leptin =    [u[7] for u in sol.u];
leptin_sen =    [u[8] for u in sol.u];
glucotoxicity = [u[9] for u in sol.u];
lipotoxicity =  [u[10] for u in sol.u];
semaglu_sub = [u[11] for u in sol.u]
semaglu_plasma = [u[12] for u in sol.u]
beta_cell = [u[13] for u in sol.u]; 
agrp = [u[14] for u in sol.u]
pomc = [u[15] for u in sol.u]
dopamine = [u[16] for u in sol.u]

if variable_name == :body_weight
    list_effect_fi_on_bw = [
        do_interpolate(effect_fi_on_bw_interp, (point / ref_food_ingestion), xs_fi_bw)
        for point in food_ingestion
    ]

    Plots.plot(sol.t, list_effect_fi_on_bw, label="effect food ingestion on body weight", xlabel= "time[days]", ylabel="effect[-]")
    Plots.plot!(twinx(),sol.t, body_weight, label="body weight", linestyle=:dash, color=:black, legend=:bottom)


elseif variable_name == :food_ingestion

    list_effect_agrp_on_fi = [
        do_interpolate(effect_agrp_on_fi_interp,(point - ref_agrp), xs_agrp_fi)
        for point in agrp
    ]        
    list_effect_pomc_on_fi = [
        do_interpolate(effect_pomc_on_fi_interp,(point - ref_pomc),xs_pomc_fi)
        for point in pomc
    ]

    list_effect_do_on_fi = [
        do_interpolate(effect_do_on_fi_interp,(point - ref_dopamine),xs_do_fi)
        for point in dopamine
    ]
    Plots.plot(sol.t, list_effect_agrp_on_fi, label="effect AgRP neuron activity on food ingestion", xlabel= "time[days]", ylabel="effect[-]")
    Plots.plot!(sol.t, list_effect_pomc_on_fi, label="effect POMC neuron activity  on food ingestion", legend=:right)
    Plots.plot!(sol.t, list_effect_do_on_fi, label="effect dopamine neuron activity  on food ingestion", legend=:right)
    Plots.plot!(twinx(),sol.t, food_ingestion, label="food ingestion", linestyle=:dash, color=:black, legend=:bottom)


elseif variable_name == :blood_glucose
    list_effect_fi_on_bg = [
        do_interpolate(effect_fi_on_bg_interp,(point / ref_food_ingestion),xs_fi_bg)
        for point in food_ingestion
    ]
    list_effect_i_on_bg = [
        do_interpolate(effect_i_on_bg_interp,(point / ref_insulin),xs_i_bg)
        for point in insulin 
    ]

    Plots.plot(sol.t, list_effect_fi_on_bg, label="effect food ingestion on blood glucose", xlabel= "time[days]", ylabel="effect[-]")
    Plots.plot!(sol.t, list_effect_i_on_bg, label="effect insulin on blood glucose")
    Plots.plot!(twinx(),sol.t, blood_glucose, label="blood glucose", linestyle=:dash, color=:black, legend=:bottom)

elseif variable_name == :glp1
    list_effect_fi_on_glp = [
        do_interpolate(effect_fi_on_glp_interp, (point / ref_food_ingestion), xs_fi_glp)
        for point in food_ingestion
    ]

    Plots.plot(sol.t, list_effect_fi_on_glp, label="effect food ingestion on GLP-1", xlabel= "time[days]", ylabel="effect[-]")
    Plots.plot!(twinx(),sol.t, glp1, label="GLP-1", linestyle=:dash, color=:black, legend=:bottom)

elseif variable_name == :insulin
    list_effect_bg_on_i = [
        do_interpolate(effect_bg_on_i_interp, (point / ref_blood_glucose), xs_bg_i)
        for point in blood_glucose
    ]
    list_effect_is_on_i = [
        do_interpolate(effect_is_on_i_interp, (point - ref_insulin_sen), xs_is_i)
        for point in insulin_sen
    ]
    
    effect_glp1_over_time = Float64[]

    for (i, bg) in enumerate(blood_glucose)
        glp = glp1[i]
        sema = semaglu_plasma[i]
        bg = blood_glucose[i]

        effect_glp_on_i1 = do_interpolate(effect_glp_on_i_interp,((glp + sema) / ref_glp_1),xs_glp_i);
        glucose_weight = clamp((bg - (90 * 1.00)) / (90 * (1.10 - 1.00)), 0.0, 1.0);
        effect_val = 1.0 + (effect_glp_on_i1-1) * glucose_weight

        push!(effect_glp1_over_time, effect_val)
    end

    list_effect_glp_on_i = [
        do_interpolate(effect_glp_on_i_interp, effect_val, xs_glp_i)
        for effect_val in effect_glp1_over_time
    ]
    list_effect_b_on_i = [
        do_interpolate(effect_b_on_i_interp, (ref_beta_cell_functioning - point), xs_b_i)

        for point in beta_cell
    ]

    Plots.plot(sol.t, list_effect_bg_on_i, label="effect blood glucose on insulin", xlabel= "time[days]", ylabel="effect[-]")
    Plots.plot!(sol.t, list_effect_is_on_i, label="effect insulin sensitivity on insulin")
    Plots.plot!(sol.t, list_effect_glp_on_i, label="effect GLP-1 + semaglutide on insulin")
    Plots.plot!(sol.t, list_effect_b_on_i, label="effect beta-cell functioning on insulin", legend =:right)

    Plots.plot!(twinx(),sol.t, insulin, label="insulin", linestyle=:dash, color=:black, legend=:topleft)

elseif variable_name == :insulin_sen
        list_effect_li_on_is = [
        do_interpolate(effect_li_on_is_interp, (point - ref_lipotoxicity), xs_li_is, ys_li_is)
        for point in lipotoxicity
    ]
    list_effect_g_on_is = [
        do_interpolate(effect_g_on_is_interp,(point - ref_glucotoxicity), xs_g_is, ys_g_is)
        for point in glucotoxicity
    ]

    Plots.plot(sol.t, list_effect_li_on_is, label="effect lipotoxicity on insulin sensitivity")
    Plots.plot!(sol.t, list_effect_g_on_is, label="effect glucotoxicity on insulin sensitivity")
    Plots.plot!(twinx(),sol.t, insulin_sen, label="insulin sensitivity", linestyle=:dash, color=:black, legend=:bottom)

elseif variable_name == :leptin
    list_effect_bw_on_l = [
        do_interpolate(effect_bw_on_l_interp, (point / ref_body_weight), xs_bw_l)
        for point in body_weight
    ]
    list_effect_ls_on_l = [
        do_interpolate(effect_ls_on_l_interp, (point - ref_leptin_sen), xs_ls_l)
        for point in leptin_sen
    ]

    Plots.plot(sol.t, list_effect_bw_on_l, label="effect body weight on leptin", xlabel= "time[days]", ylabel="effect[-]")
    Plots.plot!(sol.t, list_effect_ls_on_l, label="effect leptin sensitivity on leptin")
    Plots.plot!(twinx(),sol.t, leptin, label="leptin", linestyle=:dash, color=:black, legend=:bottom)

elseif variable_name == :leptin_sen
    list_effect_li_on_ls = [
        do_interpolate(effect_li_on_ls_interp, (point - ref_lipotoxicity), xs_li_ls)
        for point in lipotoxicity
    ]
    list_effect_l_on_ls = [
        do_interpolate(effect_l_on_ls_interp, (point / ref_leptin), xs_l_ls)
        for point in leptin
    ]

    Plots.plot(sol.t, list_effect_li_on_ls, label="effect lipotoxicity on leptin sensitivity", xlabel= "time[days]", ylabel="effect[-]")
    Plots.plot!(sol.t, list_effect_l_on_ls, label="effect leptin on leptin sensitivity")
    Plots.plot!(twinx(),sol.t, leptin_sen, label="leptin sensitivity", linestyle=:dash, color=:black, legend=:top)

elseif variable_name == :glucotoxicity
    list_effect_bg_on_g = [
        do_interpolate(effect_bg_on_g_interp, (point / ref_blood_glucose), xs_bg_g, ys_bg_g, "discrete")
        for point in blood_glucose
    ]

    Plots.plot(sol.t, list_effect_bg_on_g, label="effect blood glucose on glucotoxicity", xlabel= "time[days]", ylabel="effect[-]")
    Plots.plot!(twinx(),sol.t, glucotoxicity, label="glucotoxicity", linestyle=:dash, color=:black, legend=:topright)

elseif variable_name == :lipotoxicity
    list_effect_bw_on_li = [
        do_interpolate(effect_bw_on_li_interp, (point / ref_body_weight), xs_bw_li)
        for point in body_weight
    ]

    Plots.plot(sol.t, list_effect_bw_on_li, label="effect body weight on lipotoxicity", xlabel= "time[days]", ylabel="effect[-]")
    Plots.plot!(twinx(),sol.t, lipotoxicity, label="lipotoxicity", linestyle=:dash, color=:black, legend=:bottom)

elseif variable_name == :beta_cell
    list_effect_li_on_b = [
        do_interpolate(effect_li_on_b_interp, (point - ref_lipotoxicity), xs_li_b)
        for point in lipotoxicity
    ]
    list_effect_g_on_b = [
        do_interpolate(effect_g_on_b_interp, (point - ref_glucotoxicity), xs_g_b)
        for point in glucotoxicity
    ]

    list_effect_glp_on_b = [
        do_interpolate(effect_glp_on_b_interp, ((point+sema) / ref_glp_1), xs_glp_b)
        for (point, sema) in zip(glp1, semaglu_plasma)
    ]
    Plots.plot(sol.t, list_effect_li_on_b, label="effect lipotoxicity on beta cell functioning", xlabel= "time[days]", ylabel="effect[-]")
    Plots.plot!(sol.t, list_effect_g_on_b, label="effect glucotoxicity on beta cell functioning", xlabel= "time[days]", ylabel="effect[-]")
    Plots.plot!(sol.t, list_effect_glp_on_b, label="effect GLP-1 + semaglutide on beta cell functioning", xlabel= "time[days]", ylabel="effect[-]", legend=:top)
    Plots.plot!(twinx(),sol.t, beta_cell, label="beta cell functioning", linestyle=:dash, color=:black, legend=:bottom)

elseif variable_name == :agrp
    list_effect_i_on_agrp = [
        do_interpolate(effect_i_on_agrp_interp, (point / ref_insulin), xs_i_agrp)
        for point in insulin
    ]
    list_effect_l_on_agrp = [
        do_interpolate(effect_l_on_agrp_interp, (point / ref_leptin), xs_l_agrp)
        for point in leptin
    ]

    list_effect_do_on_agrp = [
        do_interpolate(effect_do_on_agrp_interp, (point - ref_dopamine), xs_do_agrp)
        for point in dopamine
    ]
    Plots.plot(sol.t, list_effect_i_on_agrp, label="effect insulin on AgRP neuron activity", xlabel= "time[days]", ylabel="effect[-]")
    Plots.plot!(sol.t, list_effect_l_on_agrp, label="effect leptin on AgRP neuron activity", xlabel= "time[days]", ylabel="effect[-]")
    Plots.plot!(sol.t, list_effect_do_on_agrp, label="effect dopamine on AgRP neuron activity", xlabel= "time[days]", ylabel="effect[-]",legend=:top)
    Plots.plot!(twinx(),sol.t, agrp, label="AgRP neuron activity", linestyle=:dash, color=:black, legend=:false)

elseif variable_name == :pomc
    list_effect_i_on_pomc = [
        do_interpolate(effect_i_on_pomc_interp, (point / ref_insulin), xs_i_pomc)
        for point in insulin
    ]
    list_effect_l_on_pomc = [
        do_interpolate(effect_l_on_pomc_interp, (point / ref_leptin), xs_l_pomc)
        for point in leptin
    ]

    list_effect_glp_on_pomc = [
        do_interpolate(effect_glp_on_pomc_interp, ((point+sema) / ref_glp_1), xs_glp_pomc)
        for (point, sema) in zip(glp1, semaglu_plasma)
    ]

    list_effect_do_on_pomc = [
        do_interpolate(effect_do_on_pomc_interp, (point - ref_dopamine), xs_do_pomc)
        for point in dopamine
    ]

    list_effect_agrp_on_pomc = [
        do_interpolate(effect_agrp_on_pomc_interp, (point - ref_agrp), xs_agrp_pomc)
        for point in agrp
    ]

    idx=sol.t .<= 30
    Plots.plot(sol.t[idx], list_effect_i_on_pomc[idx], label="effect insulin on POMC neuron activity", xlabel= "time[days]", ylabel="effect[-]")
    Plots.plot!(sol.t[idx], list_effect_l_on_pomc[idx], label="effect leptin on POMC neuron activity", xlabel= "time[days]", ylabel="effect[-]")
    Plots.plot!(sol.t[idx], list_effect_glp_on_pomc[idx], label="effect GLP-1 + semaglutide on POMC neuron activity", xlabel= "time[days]", ylabel="effect[-]")
    Plots.plot!(sol.t[idx], list_effect_do_on_pomc[idx], label="effect dopamine neuron activity on POMC neuron activity", xlabel= "time[days]", ylabel="effect[-]", legend=:top)
    Plots.plot!(sol.t[idx], list_effect_agrp_on_pomc[idx], label="effect AgRP neuron activity on POMC neuron activity", xlabel= "time[days]", ylabel="effect[-]")
    Plots.plot!(twinx(),sol.t[idx], pomc[idx], label="POMC neuron activity", linestyle=:dash, color=:black,legend=:false)

elseif variable_name == :dopamine
    list_effect_i_on_dopamine = [
        do_interpolate(effect_i_on_do_interp, (point / ref_insulin), xs_i_do)
        for point in insulin
    ]
    list_effect_l_on_dopamine = [
        do_interpolate(effect_l_on_do_interp, (point / ref_leptin), xs_l_do)
        for point in leptin
    ]

    list_effect_glp_on_dopamine = [
        do_interpolate(effect_glp_on_do_interp, ((point+sema) / ref_glp_1), xs_glp_do)
        for (point, sema) in zip(glp1, semaglu_plasma)
    ]

    Plots.plot(sol.t, list_effect_i_on_dopamine, label="effect insulin on dopamine neuron activity", xlabel= "time[days]", ylabel="effect[-]")
    Plots.plot!(sol.t, list_effect_l_on_dopamine, label="effect leptin on dopamine neuron activity", xlabel= "time[days]", ylabel="effect[-]")
    Plots.plot!(sol.t, list_effect_glp_on_dopamine, label="effect GLP-1 + semaglutide on dopamine neuron activity", xlabel= "time[days]", ylabel="effect[-]", legend=:top)
    Plots.plot!(twinx(),sol.t, dopamine, label="dopamine neuron activity", linestyle=:dash, color=:black, legend=false)

else
    println("Invalid variable name")
end


