%% =================================================================
% New candidate points generation (ONLY for contextual optimization)
% ==================================================================

% Build the initial value of the candidate points database chunk to be appended
if iter == 1
    db_cdpt_iter    = [sbl_seq; zeros(h_len + CDPT_F_ROWS + CDPT_G_ROWS*g_len + 2, sbl_size)];
else
    db_cdpt_iter    = [];
end

for j = 1:x_n_len
    % Calculate end points from x_n in the D-h_len cardinal directions
    % (ONLY along the axes of the INPUT variables, i.e., not contextual)
    db_cdpt_end = [repmat(x_n(:,j), 1, 2*(D-h_len));
                   zeros(CDPT_F_ROWS + CDPT_G_ROWS*g_len + 2, 2*(D-h_len))];
    n = size(db_cdpt_end, 1);
    db_cdpt_end(1:(2*n+1):end)     = 0;
%     db_cdpt_end((n+1):(2*n+1):end) = 1;

    % Append ALL samples that are not sample j, even those within the same batch of new samples
    idx_others  = true(1,X_n_len);
    idx_others(X_n_len_old+j) = false;
    db_cdpt_end = [db_cdpt_end [X_n(1:D, idx_others); zeros(CDPT_F_ROWS + CDPT_G_ROWS*g_len + 2, X_n_len-1)]];

    % Making sure that all end points are within the (normalized) bounds
    db_cdpt_end = min(1.0, max(0.0, db_cdpt_end));

    % Draw candidate points along the directions from x_n to the cdpt_ends
    for cdpt_w = (1:(B-1)) / B
        db_cdpt_iter = [db_cdpt_iter (cdpt_w*db_cdpt_end + (1-cdpt_w)*([x_n(:,j); zeros(CDPT_F_ROWS + CDPT_G_ROWS*g_len + 2, 1)]))];
    end
end

cdpts   = size(db_cdpt_iter, 2);

% Append new candidate points to the big database
db_cdpt = [db_cdpt db_cdpt_iter];
cdpts   = size(db_cdpt, 2);

if h_len
    % Project the current context information by REPLACING the corresponding rows
    db_cdpt((D-h_len+1):D,:) = repmat(h_n, 1, cdpts);
end

for cur_cdpt = 1:cdpts
    % Calculate the upper and lower bounds w.r.t. f
    % - *_ulb_lookup stores the locations of the points x_n, and the
    %   heights of the cones from x_n's
    pt_dst       = vecnorm(X_n(1:D, :) - repmat(db_cdpt(1:D, cur_cdpt), 1, X_n_len));
    f_ulb_lookup = [X_n(XN_FVAL, :); fgam*pt_dst];
    db_cdpt(CDPT_F_UB, cur_cdpt) = min(f_ulb_lookup(1, :) + f_ulb_lookup(2, :) + T_n(2, :)); % Calculating individual upper bound values TODO: should add the time-based quantities here
    db_cdpt(CDPT_F_LB, cur_cdpt) = max(f_ulb_lookup(1, :) - f_ulb_lookup(2, :) - T_n(2, :)); % Calculating individual lower bound values
    
    % Calculate the upper and lower bounds w.r.t. g
    for g_i = 1:g_len
        XN_GVAL_I = XN_GVAL + g_i - 1;
        CDPT_G_I      = CDPT_G_INFO + CDPT_G_ROWS*(g_i-1);
        g_ulb_lookup = [X_n(XN_GVAL_I, :); ggam(g_i)*pt_dst];
        db_cdpt(CDPT_G_I + CDPT_G_UB, cur_cdpt) = min(g_ulb_lookup(1, :) + g_ulb_lookup(2, :) + T_n(2+g_i, :)); % Calculating individual upper bound values
        db_cdpt(CDPT_G_I + CDPT_G_LB, cur_cdpt) = max(g_ulb_lookup(1, :) - g_ulb_lookup(2, :) - T_n(2+g_i, :)); % Calculating individual lower bound values
    end

    db_cdpt(end-1, cur_cdpt) = min(vecnorm(X_n(1:D, :) - repmat(db_cdpt(1:D, cur_cdpt), 1, X_n_len))) / diam;
end

% Calculating the uncertainties (objective and constraints) and satisfaction estimates (constraints only)
db_cdpt(CDPT_F_LAMBDA, :) = db_cdpt(CDPT_F_UB, :) - db_cdpt(CDPT_F_LB, :);
for g_i = 1:g_len
    CDPT_G_I = CDPT_G_INFO + CDPT_G_ROWS*(g_i-1);
    db_cdpt(CDPT_G_I + CDPT_G_LAMBDA, :) = db_cdpt(CDPT_G_I + CDPT_G_UB, :) - db_cdpt(CDPT_G_I + CDPT_G_LB, :);
    db_cdpt(CDPT_G_I + CDPT_G_EST, :) = ...
        sdelta*(db_cdpt(CDPT_G_I + CDPT_G_UB, :) + db_cdpt(CDPT_G_I + CDPT_G_LB, :)) / 2 + ...
        (1 - sdelta)*db_cdpt(CDPT_G_I + CDPT_G_UB, :);
end