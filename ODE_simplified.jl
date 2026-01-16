include("Data.jl")
include("Functions.jl") 

# Define interpolation functions for all effects
effect_fi_on_bw_interp = interpolate_function(xs_fi_bw, ys_fi_bw)
effect_glp_on_fi_interp = interpolate_function(xs_glp_fi, ys_glp_fi)
effect_l_on_fi_interp = interpolate_function(xs_l_fi, ys_l_fi)
effect_fi_on_bg_interp = interpolate_function(xs_fi_bg, ys_fi_bg)
effect_i_on_bg_interp = interpolate_function(xs_i_bg, ys_i_bg)
effect_fi_on_glp_interp = interpolate_function(xs_fi_glp, ys_fi_glp)
effect_bg_on_i_interp = interpolate_function(xs_bg_i, ys_bg_i_fitted)
effect_is_on_i_interp = interpolate_function(xs_is_i, ys_is_i)
effect_glp_on_i_interp = interpolate_function(xs_glp_i, ys_glp_i)
effect_i_on_is_interp = interpolate_function(xs_i_is, ys_i_is)
effect_li_on_is_interp = interpolate_function(xs_li_is, ys_li_is)
effect_g_on_is_interp = interpolate_function(xs_g_is, ys_g_is)
effect_bw_on_l_interp = interpolate_function(xs_bw_l, ys_bw_l)
effect_ls_on_l_interp = interpolate_function(xs_ls_l, ys_ls_l)
effect_li_on_ls_interp = interpolate_function(xs_li_ls, ys_li_ls)
effect_l_on_ls_interp = interpolate_function(xs_l_ls, ys_l_ls)
effect_bg_on_g_interp = interpolate_function(xs_bg_g, ys_bg_g)
effect_bw_on_li_interp = interpolate_function(xs_bw_li, ys_bw_li)
effect_b_on_i_interp = interpolate_function(xs_b_i, ys_b_i)
effect_li_on_b_interp = interpolate_function(xs_li_b, ys_li_b)
effect_g_on_b_interp = interpolate_function(xs_g_b, ys_g_b)
effect_glp_on_b_interp = interpolate_function(xs_glp_b, ys_glp_b)

function ODEs!(du,u,p,t)
     """Function contains the differential equations for the state variables.
    
    inputs:
    - du: vector with derivative of state variables
    - u: vector with current values of state variables
    - p: vector with;
        - reference values for all stock variables
        - abosrption rate constant and elimination rate constant 
    - t: vector with timepoints to simulate
    
    modifies du
    """
    # unpack p (reference values and rate constants)
    ref_body_weight, ref_food_ingestion, ref_blood_glucose,ref_glp_1, ref_insulin, ref_insulin_sen, ref_leptin, ref_leptin_sen, ref_glucotoxicity, ref_lipotoxicity, ref_semaglutide_subcutaneous, ref_beta_cell_functioning, ref_agrp, ref_pomc, ref_dopamine, absorption_rate, elimination_rate = p;

    # unpack u (current state variables)
    current_body_weight, current_food_ingestion, current_blood_glucose, current_glp_1, current_insulin, current_insulin_sen, current_leptin, current_leptin_sen, current_glucotoxicity, current_lipotoxicity, current_semaglutide_subcutaneous, current_semaglutide_plasma, current_beta_cell_functioning, current_agrp, current_pomc, current_dopamine= u;

    # Calculate variable values
    # Body weight
    target_body_weigh = target_body_weight(ref_body_weight, ref_food_ingestion, current_food_ingestion, effect_fi_on_bw_interp, xs_fi_bw)
    du[1] = (target_body_weigh - current_body_weight) / body_weight_adj_time;
    
    # Food ingestion
    target_food_ingest = target_food_ingestion_simplified(ref_food_ingestion, ref_glp_1, ref_leptin, (current_glp_1+ current_semaglutide_plasma), current_leptin, effect_glp_on_fi_interp, effect_l_on_fi_interp, xs_glp_fi, xs_l_fi)
    du[2] = (target_food_ingest - current_food_ingestion) / food_ingestion_adj_time;

    # Blood gluocose (maximal 180 mg/dL to prevent unrealistic values)
    target_blood_gluc = target_blood_glucose(ref_blood_glucose, ref_food_ingestion, ref_insulin, current_food_ingestion, current_insulin, effect_fi_on_bg_interp, effect_i_on_bg_interp, xs_fi_bg, xs_i_bg)
    target_blood_gluc = min(target_blood_gluc, 180.0);
    du[3] = (target_blood_gluc - current_blood_glucose) / blood_glucose_adj_time;

    # GLP-1 levels
    target_glp = target_glp1(ref_glp_1, ref_food_ingestion, current_food_ingestion, effect_fi_on_glp_interp, xs_fi_glp)
    du[4] = (target_glp - current_glp_1) / glp_1_adj_time;

    # Insulin levels (glucose-dependent)
    target_insul = target_insulin(ref_insulin, ref_blood_glucose, ref_insulin_sen, ref_glp_1,ref_beta_cell_functioning, current_blood_glucose,current_insulin_sen,(current_glp_1+ current_semaglutide_plasma),current_beta_cell_functioning, effect_bg_on_i_interp,effect_is_on_i_interp,effect_glp_on_i_interp,effect_b_on_i_interp, xs_bg_i, xs_is_i, xs_glp_i, xs_b_i);
    du[5] = (target_insul - current_insulin) / insulin_adj_time;
    
    # Insulin sensitivity (minimal 0.0 )
    target_insulin_sen = target_insulin_sensitivity(ref_insulin_sen, ref_lipotoxicity, ref_glucotoxicity,current_lipotoxicity,current_glucotoxicity,effect_li_on_is_interp,effect_g_on_is_interp,xs_li_is,xs_g_is)
    target_insulin_sen = max(target_insulin_sen, 0.0);
    du[6] = (target_insulin_sen - current_insulin_sen) / insulin_sen_adj_time;

    # Leptin levels
    target_lept = target_leptin(ref_leptin, ref_body_weight, ref_leptin_sen, current_body_weight, current_leptin_sen, effect_bw_on_l_interp, effect_ls_on_l_interp, xs_bw_l, xs_ls_l);
    du[7] = (target_lept - current_leptin) / leptin_adj_time;

    # Leptin sensitivity (minimal 0.0)
    target_leptin_sen = target_leptin_sensitivity(ref_leptin_sen, ref_lipotoxicity, ref_leptin, current_lipotoxicity, current_leptin, effect_li_on_ls_interp, effect_l_on_ls_interp, xs_li_ls, xs_l_ls);
    target_leptin_sen = max(target_leptin_sen, 0.0);
    du[8] = (target_leptin_sen - current_leptin_sen) / leptin_sen_adj_time;

    # Glucotoxicity (minimal 0.0)
    target_glucotox = target_glucotoxicity(ref_glucotoxicity, ref_blood_glucose, current_blood_glucose, effect_bg_on_g_interp, xs_bg_g, ys_bg_g);
    target_glucotox = max(target_glucotox, 0.0);
    du[9] = (target_glucotox - current_glucotoxicity) / glucotoxicity_adj_time;
    
    # Lipotoxicity (minimal 0.0)
    target_lipotox = target_lipotoxicity(ref_lipotoxicity, ref_body_weight, current_body_weight, effect_bw_on_li_interp, xs_bw_li)
    target_lipotox = max(target_lipotox, 0.0);
    du[10] = (target_lipotox - current_lipotoxicity) / lipotoxicity_adj_time;
    
    # Semagluide levels in plasma (11) and subcutanous (12)
    du[11] = current_semaglutide_subcutaneous * - absorption_rate;
    du[12] = current_semaglutide_subcutaneous *   absorption_rate - current_semaglutide_plasma * elimination_rate;

    # Beta-cell function
    target_beta_cell_function =target_beta_cell_functioning(ref_beta_cell_functioning, ref_lipotoxicity, ref_glucotoxicity, ref_glp_1, current_lipotoxicity, current_glucotoxicity, current_semaglutide_plasma, current_glp_1, effect_li_on_b_interp, effect_g_on_b_interp, xs_li_b, xs_g_b, xs_glp_b)
    du[13] = (target_beta_cell_function- current_beta_cell_functioning) / beta_cell_functioning_adj_time;
    
    # AgRP, POMC and dopamine neuron activity are not included
end
