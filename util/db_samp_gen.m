function [db_samp_iter, db_samp_cn_iter] = db_samp_gen(x)
    % needed information: g_len, X_n, X_n_len, x_n, x_n_len, T_n, D, fgam, ggam
    db_samp_iter    = NaN(4+4*g_len, X_n_len);
    db_samp_cn_iter = NaN(2+2*g_len, X_n_len);
    for j = 1:X_n_len
        idx_others    = true(1,X_n_len);
        idx_others(j) = false;
        pt_dst       = vecnorm(X_n(1:D, idx_others) - repmat(X_n(1:D,j), 1, X_n_len-1));
        f_ulb_lookup = [X_n(XN_FVAL, 1:end-1); fgam*pt_dst];
        [~, ub_idx] = min(f_ulb_lookup(1, :) + f_ulb_lookup(2, :) + T_n(2, 1:end-1)); % Calculating individual upper bound values
        [~, lb_idx] = max(f_ulb_lookup(1, :) - f_ulb_lookup(2, :) - T_n(2, 1:end-1)); % Calculating individual lower bound values        
        db_samp_iter(SAMP_F_UB_VTX:SAMP_F_LB_HGT) = [f_ulb_lookup(1:2, ub_idx);
                                                     f_ulb_lookup(1:2, lb_idx)];
        db_samp_cn_iter([CDPTCN_UB CDPTCN_LB]) = [X_n(XN_ID, ub_idx) X_n(XN_ID, lb_idx)]'; % Populating with the index of the cone generator samples
    
        % Calculate the upper and lower bounds w.r.t. g
        for g_i = 1:g_len
            XN_GVAL_I = XN_GVAL + g_i - 1;
            SAMP_G_I  = SAMP_G_INFO + 4*(g_i-1);
            g_ulb_lookup = [X_n(XN_GVAL_I, 1:end-1); ggam(g_i)*pt_dst];
            [~, ub_idx] = min(g_ulb_lookup(1, :) + g_ulb_lookup(2, :) + T_n(2+g_i, 1:end-1)); % Calculating individual upper bound values
            [~, lb_idx] = max(g_ulb_lookup(1, :) - g_ulb_lookup(2, :) - T_n(2+g_i, 1:end-1)); % Calculating individual lower bound values
            db_samp_iter(SAMP_G_I + (SAMP_G_UB_VTX:SAMP_G_LB_HGT)) = [g_ulb_lookup(1:2, ub_idx);
                                                                      g_ulb_lookup(1:2, lb_idx)];
            db_samp_cn_iter(2*g_i + [CDPTCN_UB CDPTCN_LB]) = [X_n(XN_ID, ub_idx) X_n(XN_ID, lb_idx)]'; % Populating with the index of the cone generator samples
        end
    end