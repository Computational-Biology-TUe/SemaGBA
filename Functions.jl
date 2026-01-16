using Interpolations

function run_simulation(inputfile, odefile; callback_type= "none")
    """ 
    Run an ODE simulation for a given input scenario
    
    Inputs:
    - inputfile: name of julia file defining;
        - Initial and reference values for all variables
        - Adjustment time of all variables
        - Absorption rate constant and elimination rate constant
        - Duration of simulation 
        - Injection doses corrected with bioavailability and distribution volume
        - Injection timepoints  
    - odefile: name of a julia file defining the ODE system
    - callback: type of callback that should be applied
        - 'none': no callbacks
        - 'semaglutide': callback to simulate semaglutide treatment. 
                        Dose and timepoints are specified in inputfile
        - 'overeating': callback to simulate overeating
        - 'overeating with semaglutide': callback to simulate semaglutide treatment during overeating

    Returns:
    - sol: ODESolution object containing the time vector and state variables
    """
    include(inputfile)
    include(odefile)
    
    # Define ODE
    tspan = (0.0, t_days)
    u0 = [initial_body_weight, initial_food_ingestion, initial_blood_glucose, 
        initial_glp_1,  initial_insulin, initial_insulin_sen, 
        initial_leptin, initial_leptin_sen, initial_glucotoxicity, 
        initial_lipotoxicity, initial_semaglutide_subcutaneous, initial_semaglutide_plasma, 
        initial_beta_cell_functioning, initial_agrp, initial_pomc, initial_dopamine];

    p = (ref_body_weight, ref_food_ingestion, ref_blood_glucose,
        ref_glp_1, ref_insulin, ref_insulin_sen, 
        ref_leptin, ref_leptin_sen, ref_glucotoxicity,
        ref_lipotoxicity, ref_semaglutide_subcutaneous, 
        ref_beta_cell_functioning, 
        ref_agrp, ref_pomc,  ref_dopamine,
        absorption_rate, elimination_rate);
    
    # Concatenate injection arrays 
    injections_time_points = vcat(release_times...)
    injections_doses = vcat([fill(dose, length(t)) for (dose,t) in zip(injection_doses, release_times)]...)
    
    # Determine callback
    if callback_type == "none"
        cb = CallbackSet()

    elseif callback_type == "semaglutide"        
        cb = CallbackSet(semaglutide_callback(injections_time_points, injections_doses))

    elseif callback_type == "overeating"
        cb = CallbackSet(overeating_callback())

    elseif callback_type == "overeating with semaglutide"
        cb_sema = semaglutide_callback(injections_time_points, injections_doses)
        cb_food = overeating_callback()
        cb = CallbackSet(cb_food,cb_sema);
    else 
        error("Unknown type. Use \"none\",\"semaglutide\", \"overeating\"  or \"overeating with semaglutide\".")
    end 

    # Run simulation
    prob = ODEProblem(ODEs!, u0, tspan, p)
    sol = solve(prob, RadauIIA5(); tstops=injections_time_points, saveat=0.1, callback=cb)

    return sol
end

"""
Callback functions to add semaglutide injections as external input
"""

function semaglutide_callback(injections_time_points, injections_doses)
    condition = (u,t,integrator) -> any(isapprox(t, timepoint; atol=1e-8, rtol=0) for timepoint in injections_time_points)
    function add_semaglutide!(integrator)
        t = integrator.t
        idx = findfirst(timepoint -> isapprox(t, timepoint; atol=1e-8, rtol=0), injections_time_points)
        if idx !== nothing && 1 <= idx <= length(injections_doses)
            integrator.u[11] += injections_doses[idx]
        else
            @warn "No matching injection or index out of range for t=$t" idx
        end
    end
    return DiscreteCallback(condition, add_semaglutide!)
end

"""
Callback functions to keep food ingestion constant during the simulation
"""
function overeating_callback()
    # Apply the callback during the period that food ingestion should remain constant
    condition = (u,t,integrator) -> t < 365*5
    function keep_overeating!(integrator)
        # Set food intake to a constant value
        integrator.u[2] = 3000
    end
    return DiscreteCallback(condition, keep_overeating!)
end


"""
Functions to calculate used to calculate target values
"""

function interpolate_function(x_data,y_data)
    """
    Creates an interpolation function from the given x and y data points
    
    Inputs:
    - x_data : array of x-values.
    - y_data : array of y-values, corresponding to 'x_data'.

    Returns a callable function that estimates intermediate y-values based on linear interpolation.
    """
    return linear_interpolation(x_data,y_data)
end

function do_interpolate(interpolate_function, point, x_data, y_data=[], type = "continuous")
    """ 
    Returns the y-value corresponding to a given x-value ('point'), using continuous interpolation or discrete lookup.
    
    Inputs:
    - interpolate_function: a function created by interpolate_function(...), used for continuous interpolation.
    - point: the x-value for which the corresponding y-value is requested.
    - x_data: array of x-values.
    - y_data: array of y-values, corresponding to 'x_data'. Only required for type 'discrete'.
    - type: Method for value determination.
        - "continuous": uses interpolation to compute the intermediate value. 'Point' is clamped to the range of 'x_data'
        - "discrete": returns the y-value corresponding to the x-value that is closest to 'point'.
    
    Returns a single interpolated or nearest y-value
    """    
    point = float(point)

    # Determine minimum and maximum of the x-values
    x_min = minimum(x_data)
    x_max = maximum(x_data)

    if type == "continuous"
        # Clamp point to valid range to prevent extrapolation
        if point < x_min
            point = x_min
        elseif point > x_max
            point = x_max
        end
        # Use interpolation function to get y-value
        return interpolate_function(point)
    
    elseif type == "discrete"
        # Check if y_data is provided
        if y_data == []
            error("list of y-values is required")
        else
        # Find index of x-value closest to point 
        idx_closest = argmin(abs.(x_data .- point)) 
        # Return corresponding y-value
        return y_data[idx_closest]
        end
    else
        error("Unknown type. Use \"continuous\" or \"discrete\".")
    end
end

function calculate_target_variable(reference_value, effects, type = "product")
    """
    Calculates the target value for given variable 
    by multiplying the effects of other variables on this variable, 
    and the reference value of this variable.

    Inputs:
    - reference_value: constant
    - effects: value that indicates the effect of one variable on another at a timepoint
    - type; 'product', 'sum', or 'sum product'; determines how to combine the effects
        - 'product': multiply all effects and reference value
        - 'sum': sum all effects and reference value
        - 'sum product ': average the effects and multiply by reference value
    returns the product of all effects for the given variable
    """
    if type == "product"
        return prod(e for e in effects) * reference_value
    elseif type == "sum"
        return (sum(e for e in effects)) + reference_value
    elseif type == "sum product"
       n = length(effects)
        if n == 0
            return reference_value
        else
            return (sum(e for e in effects) /n) * reference_value
        end
    else
        error("unknown type. Use 'product', 'sum' or 'sum product'.")
    end
end


"""
Functions to calculate the target values for each the state variables
"""

 # Body weigt, influenced by food ingestion
function target_body_weight(ref_body_weight, ref_food_ingestion, current_food_ingestion, effect_fi_on_bw_interp, xs_fi_bw)
    effect_fi_on_bw = do_interpolate(effect_fi_on_bw_interp,(current_food_ingestion / ref_food_ingestion), xs_fi_bw);
    target_bw = calculate_target_variable(ref_body_weight, effect_fi_on_bw);
end

# Food ingestion for the simplified model, influenced by GLP-1 and leptin levels
function target_food_ingestion_simplified(ref_food_ingestion, ref_glp_1, ref_leptin, current_glp_1, current_leptin, effect_glp_on_fi_interp, effect_l_on_fi_interp, xs_glp_fi, xs_l_fi)
    effect_glp_on_fi = do_interpolate(effect_glp_on_fi_interp,(current_glp_1 / ref_glp_1),xs_glp_fi);
    effect_l_on_fi = do_interpolate(effect_l_on_fi_interp,(current_leptin / ref_leptin),xs_l_fi);
    
    target_fi = calculate_target_variable(ref_food_ingestion, [effect_glp_on_fi, effect_l_on_fi]);
end

# Food ingestion for the neural model, influenced by AgRP, POMC, and dopamine neural activity
function target_food_ingestion_neural(ref_food_ingestion, ref_agrp, ref_pomc, ref_dopamine, current_agrp, current_pomc, current_dopamine, effect_agrp_on_fi_interp, effect_pomc_on_fi_interp, effect_do_on_fi_interp, xs_agrp_fi, xs_pomc_fi, xs_do_fi)
    effect_agrp_on_fi = do_interpolate(effect_agrp_on_fi_interp,(current_agrp - ref_agrp),xs_agrp_fi);
    effect_pomc_on_fi = do_interpolate(effect_pomc_on_fi_interp,(current_pomc - ref_pomc),xs_pomc_fi);
    effect_do_on_fi = do_interpolate(effect_do_on_fi_interp,(current_dopamine - ref_dopamine),xs_do_fi);

    target_fi = calculate_target_variable(ref_food_ingestion, [effect_agrp_on_fi, effect_pomc_on_fi, effect_do_on_fi]);
end


# Blood glucose, influenced by food ingestion, and insulin
function target_blood_glucose(ref_blood_glucose, ref_food_ingestion, ref_insulin, current_food_ingestion, current_insulin, effect_fi_on_bg_interp, effect_i_on_bg_interp, xs_fi_bg, xs_i_bg)
    effect_fi_on_bg = do_interpolate(effect_fi_on_bg_interp,(current_food_ingestion / ref_food_ingestion),xs_fi_bg);
    effect_i_on_bg = do_interpolate(effect_i_on_bg_interp,(current_insulin / ref_insulin),xs_i_bg);
    
    target_bg = calculate_target_variable(ref_blood_glucose, [effect_fi_on_bg, effect_i_on_bg]);
end

# GLP-1, influenced by food ingestion
function target_glp1(ref_glp_1, ref_food_ingestion, current_food_ingestion, effect_fi_on_glp_interp, xs_fi_glp);
    effect_fi_on_glp = do_interpolate(effect_fi_on_glp_interp,(current_food_ingestion / ref_food_ingestion),xs_fi_glp);

    target_glp1 = calculate_target_variable(ref_glp_1, [effect_fi_on_glp]);
end


# Insulin, influenced by blood glucucose, insulin sensitivity, beta-cell function, and GLP-1 (dependent on blood glucose), 
function target_insulin(ref_insulin, ref_blood_glucose, ref_insulin_sen, ref_glp_1,ref_beta_cell_functioning, current_blood_glucose,current_insulin_sen,current_glp_1,current_beta_cell_functioning, effect_bg_on_i_interp,effect_is_on_i_interp,effect_glp_on_i_interp,effect_b_on_i_interp, xs_bg_i, xs_is_i, xs_glp_i, xs_b_i);
    effect_bg_on_i = do_interpolate(effect_bg_on_i_interp,(current_blood_glucose/ ref_blood_glucose),xs_bg_i);
    effect_is_on_i = do_interpolate(effect_is_on_i_interp,(current_insulin_sen - ref_insulin_sen),xs_is_i);
    effect_b_on_i = do_interpolate(effect_b_on_i_interp,(ref_beta_cell_functioning - current_beta_cell_functioning),xs_b_i);

    # The effect of GLP-1 is dependent on blood glucose levels
    if current_blood_glucose == ref_blood_glucose
        effect_glp1 = 1
    elseif (current_blood_glucose / 90 > 1.05 ) && (current_glp_1/ref_glp_1 > 1.01);
        effect_glp1= ((current_glp_1 / ref_glp_1)* (current_blood_glucose/90));
    else 
        effect_glp1 = 1
    end
    effect_glp_on_i = do_interpolate(effect_glp_on_i_interp,effect_glp1,xs_glp_i);

    target_insulin = calculate_target_variable(ref_insulin, [effect_bg_on_i, effect_glp_on_i, effect_is_on_i, effect_b_on_i]);
end

# Insulin sensitivity, influenced by lipotoxicity and glucotoxicity
function target_insulin_sensitivity(ref_insulin_sen, ref_lipotoxicity, ref_glucotoxicity,current_lipotoxicity,current_glucotoxicity,effect_li_on_is_interp,effect_g_on_is_interp,xs_li_is,xs_g_is)
    effect_li_on_is = do_interpolate(effect_li_on_is_interp,(current_lipotoxicity - ref_lipotoxicity),xs_li_is, ys_li_is);
    effect_g_on_is = do_interpolate(effect_g_on_is_interp,(current_glucotoxicity - ref_glucotoxicity),xs_g_is, ys_g_is);

    target_insulin_sen = calculate_target_variable(ref_insulin_sen, [effect_li_on_is, effect_g_on_is] ,"sum product");
end

# Leptin, influenced by body weight and leptin sensitivity
function target_leptin(ref_leptin, ref_body_weight, ref_leptin_sen, current_body_weight, current_leptin_sen, effect_bw_on_l_interp, effect_ls_on_l_interp, xs_bw_l, xs_ls_l);
    
    effect_bw_on_l = do_interpolate(effect_bw_on_l_interp,(current_body_weight / ref_body_weight),xs_bw_l);
    effect_ls_on_l = do_interpolate(effect_ls_on_l_interp,(current_leptin_sen - ref_leptin_sen),xs_ls_l);
    
    target_leptin = calculate_target_variable(ref_leptin, [effect_bw_on_l, effect_ls_on_l]);
end

# Leptin sensitivity, influenced by lipotoxicity and leptin
function target_leptin_sensitivity(ref_leptin_sen, ref_lipotoxicity, ref_leptin, current_lipotoxicity, current_leptin, effect_li_on_ls_interp, effect_l_on_ls_interp, xs_li_ls, xs_l_ls)
    effect_li_on_ls = do_interpolate(effect_li_on_ls_interp,(current_lipotoxicity - ref_lipotoxicity),xs_li_ls);
    effect_l_on_ls = do_interpolate(effect_l_on_ls_interp,(current_leptin / ref_leptin),xs_l_ls);
    
    target_leptin_sen = calculate_target_variable(ref_leptin_sen, [effect_li_on_ls, effect_l_on_ls] ,"sum product");
end

# Glucotoxicity, influenced by blood glucose 
function target_glucotoxicity(ref_glucotoxicity, ref_blood_glucose, current_blood_glucose, effect_bg_on_g_interp, xs_bg_g, ys_bg_g);   
   effect_bg_on_g = do_interpolate(effect_bg_on_g_interp,(current_blood_glucose / ref_blood_glucose),xs_bg_g, ys_bg_g ,"discrete");

   target_glucotox = calculate_target_variable(ref_glucotoxicity, [effect_bg_on_g],"sum");
end

# Lipotoxicity, influenced by body weight
function target_lipotoxicity(ref_lipotoxicity, ref_body_weight, current_body_weight, effect_bw_on_li_interp, xs_bw_li);
    effect_bw_on_li = do_interpolate(effect_bw_on_li_interp,(current_body_weight / ref_body_weight),xs_bw_li, ys_bw_li ,"discrete");

    target_lipotox = calculate_target_variable(ref_lipotoxicity, [effect_bw_on_li] ,"sum");
end

# Beta-cell functioning, incluenced by lipotoxicity, glucotoxicity and GLP-1
function target_beta_cell_functioning(ref_beta_cell_functioning, ref_lipotoxicity, ref_glucotoxicity, ref_glp_1, current_lipotoxicity, current_glucotoxicity, current_semaglutide_plasma, current_glp_1, effect_li_on_b_interp, effect_g_on_b_interp, xs_li_b, xs_g_b, xs_glp_b);
    effect_li_b = do_interpolate(effect_li_on_b_interp, (current_lipotoxicity - ref_lipotoxicity), xs_li_b);
    effect_g_b = do_interpolate(effect_g_on_b_interp, (current_glucotoxicity - ref_glucotoxicity), xs_g_b);
    effect_glp_b = do_interpolate(effect_glp_on_b_interp,((current_semaglutide_plasma +  current_glp_1) / ref_glp_1), xs_glp_b);
    
    target_beta = calculate_target_variable(ref_beta_cell_functioning, [effect_li_b, effect_g_b, effect_glp_b] ,"sum product");
end

# AgRP neuron activity, influenced by insulin, leptin, and dopamine neuron activity
function target_AgRP(ref_agrp, ref_dopamine, ref_insulin, ref_leptin, current_dopamine, current_insulin, current_leptin, effect_do_on_agrp_interp, effect_i_on_agrp_interp, effect_l_on_agrp_interp, xs_do_agrp, xs_i_agrp, xs_l_agrp)
    effect_i_agrp = do_interpolate(effect_i_on_agrp_interp, (current_insulin / ref_insulin), xs_i_agrp);
    effect_l_agrp = do_interpolate(effect_l_on_agrp_interp, (current_leptin / ref_leptin), xs_l_agrp);
    effect_do_agrp = do_interpolate(effect_do_on_agrp_interp, (current_dopamine - ref_dopamine), xs_do_agrp);

    target_agrp = calculate_target_variable(ref_agrp, [effect_i_agrp, effect_l_agrp, effect_do_agrp], "sum product");
end

# POMC neuron activity, influenced by insulin, leptin, GLP-1, dopamine neuron activity and AgRP neuron activity
function target_POMC(ref_pomc, ref_glp_1, ref_dopamine, ref_insulin, ref_leptin, ref_agrp, current_glp_1, current_dopamine, current_insulin, current_leptin, current_agrp, effect_i_on_pomc_interp, effect_l_on_pomc_interp, effect_glp_on_pomc_interp, effect_do_on_pomc_interp, effect_agrp_on_pomc_interp, xs_i_pomc, xs_l_pomc, xs_glp_pomc, xs_do_pomc, xs_agrp_pomc);
    effect_i_pomc = do_interpolate(effect_i_on_pomc_interp, (current_insulin / ref_insulin), xs_i_pomc);
    effect_l_pomc = do_interpolate(effect_l_on_pomc_interp, (current_leptin / ref_leptin), xs_l_pomc);
    effect_glp_pomc = do_interpolate(effect_glp_on_pomc_interp, (current_glp_1 / ref_glp_1), xs_glp_pomc);
    effect_do_pomc = do_interpolate(effect_do_on_pomc_interp, (current_dopamine - ref_dopamine), xs_do_pomc);
    effect_agrp_pomc = do_interpolate(effect_agrp_on_pomc_interp, (current_agrp - ref_agrp), xs_agrp_pomc);

    target_pomc = calculate_target_variable(ref_pomc, [effect_glp_pomc, effect_do_pomc, effect_i_pomc, effect_l_pomc, effect_agrp_pomc], "sum product");
end

# Doapmine neuron activity, influenced by insuline, leptin and GLP-1
function target_Dopamine(ref_dopamine, ref_insulin, ref_leptin, ref_glp_1, current_insulin, current_leptin, current_glp_1, effect_i_do_interp, effect_l_on_do_interp, effect_glp_do_interp, xs_i_do, xs_l_do, xs_glp_do)
    effect_i_dopamine = do_interpolate(effect_i_on_do_interp, (current_insulin / ref_insulin), xs_i_do);
    effect_l_dopamine = do_interpolate(effect_l_on_do_interp, (current_leptin / ref_leptin), xs_l_do);
    effect_glp_dopamine = do_interpolate(effect_glp_on_do_interp, (current_glp_1 / ref_glp_1), xs_glp_do);

    target_dopamine = calculate_target_variable(ref_dopamine, [effect_i_dopamine, effect_l_dopamine, effect_glp_dopamine], "sum product");
end

"""
Functions for the barplot
"""

function normalized_change(inputfile, odefile; callback_type = "none")
    include(inputfile)

    sol = run_simulation(inputfile, odefile, callback_type = callback_type)                 

    body_weight     = ([round(u[1], digits=12) for u in sol.u][end] -ref_body_weight)/ref_body_weight 
    food_ingestion   = ([round(u[2], digits=12) for u in sol.u][end] -ref_food_ingestion)/ref_food_ingestion 
    blood_glucose   = ([round(u[3], digits=12) for u in sol.u][end] - ref_blood_glucose) / ref_blood_glucose 
    glp1            = ([round(u[4], digits=12) for u in sol.u][end]  - ref_glp_1) / ref_glp_1 
    insulin          = ([round(u[5], digits=12) for u in sol.u][end]  - ref_insulin) / ref_insulin 
    leptin           = ([round(u[7], digits=12) for u in sol.u][end]  - ref_leptin) / ref_leptin 
    insulin_sen      = ([round(u[6], digits=12) for u in sol.u][end]  - ref_insulin_sen)
    leptin_sen       = ([round(u[8], digits=12) for u in sol.u][end]  - ref_leptin_sen) 
    beta             = ([round(u[13], digits=12) for u in sol.u][end]  - ref_beta_cell_functioning) 
    glucotoxicity    = ([round(u[9], digits=12) for u in sol.u][end]  - ref_glucotoxicity) 
    lipotoxicity     = ([round(u[10], digits=12) for u in sol.u][end]  - ref_lipotoxicity) 
    
    return [body_weight, food_ingestion, blood_glucose, glp1, insulin, leptin, insulin_sen, leptin_sen, beta, glucotoxicity, lipotoxicity]
end


function compute_x(x, dodge; width=1, gap=0.2, dodge_gap=0.03)
    """
    Calculates x-positions of grouped bars in bar plot

    Inputs: 
    - x = original position of all variables
    - dodge = category of the bar
    - width = width of the bar
    - gap = space between groups
    - dodge_gap = space between bars of different categories
    
    Returns:
    - shifted position of all variables
    """
    
    scale_width(dodge_gap, n_dodge) = (1 - (n_dodge - 1) * dodge_gap) / n_dodge
    function shift_dodge(i, dodge_width, dodge_gap)
        (dodge_width - 1) / 2 + (i - 1) * (dodge_width + dodge_gap)
    end
    width *= 1 - gap
    n_dodge = maximum(dodge)
    dodge_width = scale_width(dodge_gap, n_dodge)
    shifts = shift_dodge.(dodge, dodge_width, dodge_gap)
    return x .+ width .* shifts
end

