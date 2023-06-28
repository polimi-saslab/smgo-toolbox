% ===================================
% For candidate points database
% ===================================
db_cdpt([CDPT_F_UB_HGT CDPT_F_LB_HGT],:) = ...
    (fgam / fgam_prev) .* db_cdpt([CDPT_F_UB_HGT CDPT_F_LB_HGT],:);        
for g_i = 1:g_len
    CDPT_G_I = CDPT_G_INFO + CDPT_G_ROWS*(g_i-1);
    db_cdpt(CDPT_G_I + [CDPT_G_UB_HGT CDPT_G_LB_HGT],:) = ...
        (ggam(g_i) / ggam_prev(g_i)) .* db_cdpt(CDPT_G_I + [CDPT_G_UB_HGT CDPT_G_LB_HGT],:);
end

cdpts = size(db_cdpt, 2);

% T_n_coll_ub (_lb) is a TEMPORARY [1 x n] matrix
T_n_coll_ub  = NaN(1, size(db_cdpt,2));
T_n_coll_lb  = NaN(1, size(db_cdpt,2));
for i = 1:X_n_len
    % X_n(XN_ID,i) is the sample ID of sample i
    % For each sample i, the next two lines update the "stacked uncertainties"
    % for the upper (lower) bound values for those samples for which sample i
    % generated their upper (lower) bound. Please also see documentation to 
    % T_n (smgo.m) and db_samp_cn (update_gen_cdpt.m)
    T_n_coll_ub(db_cdpt_cn(CDPTCN_UB,:)==X_n(XN_ID,i)) = T_n(2,i);
    T_n_coll_lb(db_cdpt_cn(CDPTCN_LB,:)==X_n(XN_ID,i)) = T_n(2,i);
end

for j = 1:x_n_len
    xn_id_tmp = X_n(XN_ID, X_n_len_old+j);
    z_tmp     = X_n(XN_FVAL, X_n_len_old+j);

    % Calculating if I should update the f upper (lower) bounds
    new_cone_hgt = fgam * vecnorm(db_cdpt(1:D,:) - repmat(x_n(:,j), 1, cdpts));
    new_f_ub     = z_tmp + feps + new_cone_hgt; 
    new_f_lb     = z_tmp - feps - new_cone_hgt;
    
    % f_ub_update (_lb_) is a 'bool' vector selecting the columns in which the new upper (lower) bound is TIGHTER than the existing one
    f_ub_update  = (db_cdpt(CDPT_F_UB_VTX,:) + db_cdpt(CDPT_F_UB_HGT,:) + feps + T_n_coll_ub) > new_f_ub;
    f_lb_update  = (db_cdpt(CDPT_F_LB_VTX,:) - db_cdpt(CDPT_F_LB_HGT,:) - feps - T_n_coll_lb) < new_f_lb;
    
    % Replacing the f upper (lower) bounds data if needed
    db_cdpt(CDPT_F_UB_VTX, f_ub_update) = z_tmp;
    db_cdpt(CDPT_F_UB_HGT, f_ub_update) = new_cone_hgt(:, f_ub_update);
    db_cdpt(CDPT_F_LB_VTX, f_lb_update) = z_tmp;
    db_cdpt(CDPT_F_LB_HGT, f_lb_update) = new_cone_hgt(:, f_lb_update);
    db_cdpt_cn(CDPTCN_UB, f_ub_update) = xn_id_tmp; % Update the ID of the new cone-generator sample
    db_cdpt_cn(CDPTCN_LB, f_lb_update) = xn_id_tmp;

    for g_i = 1:g_len
        CDPT_G_I = CDPT_G_INFO + CDPT_G_ROWS*(g_i-1);
        c_tmp    = X_n(XN_GVAL+g_i-1, X_n_len_old+j);

        % Calculating if I should update the g upper (lower) bounds
        new_cone_hgt = ggam(g_i) * vecnorm(db_cdpt(1:D,:) - repmat(x_n(:,j), 1, cdpts));
        new_g_ub     = c_tmp + geps(g_i) + new_cone_hgt; 
        new_g_lb     = c_tmp - geps(g_i) - new_cone_hgt;
    
        T_n_coll_ub  = NaN(1, size(db_cdpt,2));
        T_n_coll_lb  = NaN(1, size(db_cdpt,2));
        for i = 1:X_n_len
            T_n_coll_ub(db_cdpt_cn(2*g_i+CDPTCN_UB,:)==X_n(XN_ID,i)) = T_n(2+g_i,i);
            T_n_coll_lb(db_cdpt_cn(2*g_i+CDPTCN_LB,:)==X_n(XN_ID,i)) = T_n(2+g_i,i);
        end
    
        g_ub_update  = new_g_ub < (db_cdpt(CDPT_G_I + CDPT_G_UB_VTX,:) + db_cdpt(CDPT_G_I + CDPT_G_UB_HGT,:) + geps(g_i) + T_n_coll_ub);
        g_lb_update  = new_g_lb > (db_cdpt(CDPT_G_I + CDPT_G_LB_VTX,:) - db_cdpt(CDPT_G_I + CDPT_G_LB_HGT,:) - geps(g_i) - T_n_coll_lb);
    
        % Replacing the g upper (lower) bounds data if needed
        db_cdpt(CDPT_G_I + CDPT_G_UB_VTX, g_ub_update) = c_tmp;
        db_cdpt(CDPT_G_I + CDPT_G_UB_HGT, g_ub_update) = new_cone_hgt(:, g_ub_update);
        db_cdpt(CDPT_G_I + CDPT_G_LB_VTX, g_lb_update) = c_tmp;
        db_cdpt(CDPT_G_I + CDPT_G_LB_HGT, g_lb_update) = new_cone_hgt(:, g_lb_update);
    
        db_cdpt_cn(2*g_i+CDPTCN_UB, g_ub_update) = xn_id_tmp; % Update the ID of the new cone-generator sample
        db_cdpt_cn(2*g_i+CDPTCN_LB, g_lb_update) = xn_id_tmp;
    end

    db_cdpt(end-1,:) = min(db_cdpt(end-1,:), vecnorm(db_cdpt(1:D,:) - repmat(x_n(:,j), 1, cdpts)) / diam);
end

% Incrementing age of existing candidate points
db_cdpt(end,:) = db_cdpt(end,:) + 1;