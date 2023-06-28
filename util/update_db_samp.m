% ======================================================
% For samples database (used only for time-varying case)
% ======================================================

% Scale the cone heights according to the change in Lipschitz constants
db_samp([SAMP_F_UB_HGT SAMP_F_LB_HGT], :) = ...
    (fgam / fgam_prev) .* db_samp([SAMP_F_UB_HGT SAMP_F_LB_HGT], :);
for g_i = 1:g_len
    SAMP_G_I = SAMP_G_INFO + 4*(g_i-1);
    db_samp(SAMP_G_I + [SAMP_G_UB_HGT SAMP_G_LB_HGT], :) = ...
        (ggam(g_i) / ggam_prev(g_i)) .* db_samp(SAMP_G_I + [SAMP_G_UB_HGT SAMP_G_LB_HGT], :);
end

for j = 1:x_n_len
    xn_id_tmp = X_n(XN_ID, X_n_len_old+j);
    z_tmp     = X_n(XN_FVAL, X_n_len_old+j);

    % Calculating if I should update the f upper (lower) bounds
    new_cone_hgt = fgam * vecnorm(X_n(1:D, 1:X_n_len_old) - repmat(x_n(:,j), 1, X_n_len_old));
    new_f_ub     = z_tmp + feps + new_cone_hgt; 
    new_f_lb     = z_tmp - feps - new_cone_hgt;

    % T_n_coll_ub (_lb) is a TEMPORARY [1 x n] matrix
    T_n_coll_ub  = NaN(1, X_n_len_old);
    T_n_coll_lb  = NaN(1, X_n_len_old);
    for i = 1:X_n_len_old
        % X_n(XN_ID,i) is the sample ID of sample i
        % For each sample i, the next two lines updates the "stacked uncertainties"
        % for the upper (lower) bound values for those samples for which sample i
        % generated their upper (lower) bound. Please also see documentation to 
        % T_n (smgo.m) and db_samp_cn (update_gen_cdpt.m)
        T_n_coll_ub(db_samp_cn(CDPTCN_UB,:)==X_n(XN_ID,i)) = T_n(2,i);
        T_n_coll_lb(db_samp_cn(CDPTCN_LB,:)==X_n(XN_ID,i)) = T_n(2,i);
    end

    % f_ub_update (_lb_) contains Bool on whether the upper (lower) bound from incoming sample is TIGHTER than the existing one.
    % The expressions to the right of < (>) are the old bounds: cone vertex value + (-) cone height + (-) stacked uncertainty
    f_ub_update = (new_f_ub < (db_samp(SAMP_F_UB_VTX, :) + db_samp(SAMP_F_UB_HGT, :) + feps + T_n_coll_ub));
    f_lb_update = (new_f_lb > (db_samp(SAMP_F_LB_VTX, :) - db_samp(SAMP_F_LB_HGT, :) - feps - T_n_coll_lb));
    
    % Replacing the f upper (lower) bounds data if marked in f_ub_update (_lb_)
    db_samp(SAMP_F_UB_VTX, f_ub_update) = z_tmp;
    db_samp(SAMP_F_UB_HGT, f_ub_update) = new_cone_hgt(:, f_ub_update);
    db_samp(SAMP_F_LB_VTX, f_lb_update) = z_tmp;
    db_samp(SAMP_F_LB_HGT, f_lb_update) = new_cone_hgt(:, f_lb_update);
    % New sample index/ID becomes the cone generator
    db_samp_cn(CDPTCN_UB, f_ub_update) = xn_id_tmp;
    db_samp_cn(CDPTCN_LB, f_lb_update) = xn_id_tmp;

    for g_i = 1:g_len
        SAMP_G_I = SAMP_G_INFO + 4*(g_i-1);
        c_tmp    = X_n(XN_GVAL+g_i-1, X_n_len_old+j);

        % Calculating if I should update the g upper (lower) bounds
        new_cone_hgt = ggam(g_i) * vecnorm(db_samp(1:D, :) - repmat(x_n(:,j), 1, size(db_samp,2)));
        new_g_ub     = c_tmp + geps(g_i) + new_cone_hgt; 
        new_g_lb     = c_tmp - geps(g_i) - new_cone_hgt;
    
        % T_n_coll_ub (_lb) is a TEMPORARY [1 x n] matrix
        T_n_coll_ub  = NaN(1, size(db_samp_cn,2));
        T_n_coll_lb  = NaN(1, size(db_samp_cn,2));
        for i = 1:X_n_len_old
            % X_n(XN_ID,i) is the sample ID of sample i
            % For each sample i, the next two lines updates the "stacked uncertainties"
            % for the upper (lower) bound values for those samples for which sample i
            % generated their upper (lower) bound. Please also see documentation to 
            % T_n (smgo.m) and db_samp_cn (update_gen_cdpt.m)
            T_n_coll_ub(db_samp_cn(2*g_i+CDPTCN_UB,:)==X_n(XN_ID,i)) = T_n(2+g_i,i);
            T_n_coll_lb(db_samp_cn(2*g_i+CDPTCN_LB,:)==X_n(XN_ID,i)) = T_n(2+g_i,i);
        end
    
        % g_ub_update (_lb_) contains Bool on whether the upper (lower) bound from incoming sample is TIGHTER than the existing one.
        % The expressions to the right of < (>) are the old bounds: cone vertex value + (-) cone height + (-) stacked uncertainty
        g_ub_update = (new_g_ub < (db_samp(SAMP_G_I + SAMP_G_UB_VTX, :) + db_samp(SAMP_G_I + SAMP_G_UB_HGT, :) + geps(g_i) + T_n_coll_ub));
        g_lb_update = (new_g_lb > (db_samp(SAMP_G_I + SAMP_G_LB_VTX, :) - db_samp(SAMP_G_I + SAMP_G_LB_HGT, :) - geps(g_i) - T_n_coll_lb));
    
        % Replacing the f upper (lower) bounds data if marked in f_ub_update (_lb_)
        db_samp(SAMP_G_I + SAMP_G_UB_VTX, g_ub_update) = c_tmp;
        db_samp(SAMP_G_I + SAMP_G_UB_HGT, g_ub_update) = new_cone_hgt(:, g_ub_update);
        db_samp(SAMP_G_I + SAMP_G_LB_VTX, g_lb_update) = c_tmp;
        db_samp(SAMP_G_I + SAMP_G_LB_HGT, g_lb_update) = new_cone_hgt(:, g_lb_update);
        % New sample index/ID becomes the cone generator
        db_samp_cn(2*g_i+CDPTCN_UB, g_ub_update) = xn_id_tmp;
        db_samp_cn(2*g_i+CDPTCN_LB, g_lb_update) = xn_id_tmp;
    end
end