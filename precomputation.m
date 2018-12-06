function [] = precomputation()
    if ~ isfile('+autogen_func/swing_model.m')
        compute_swing_model();
    end

    if ~isfile('+autogen_func/impact_model.m')
        compute_impact_model();
    end
end