%% =================================================================
% New candidate points generation (ONLY for contextual optimization)
% ==================================================================

% Build the initial value of the candidate points database chunk to be appended
if isempty(db_cdpt)
    db_cdpt_iter    = [zeros(h_len,sbl_size); sbl_seq; zeros(CDPT_F_ROWS + CDPT_G_ROWS*g_len + 2, sbl_size)];
    for j = 1:X0_len
        % Calculate end points in the cardinal directions
        db_cdpt_end = [repmat(X0(1:D,j), 1, 2*D);
                       zeros(CDPT_F_ROWS + CDPT_G_ROWS*g_len + 2, 2*D)];
        n = size(db_cdpt_end, 1);
        db_cdpt_end(1:(2*n+1):end)     = 0;
        db_cdpt_end((n+1):(2*n+1):end) = 1;
    
        % Check if each end point generated along the cardinal direction is within the inequality constraints
        % If outside these constraints, project to the boundary
        if isfield(options,'ineq')
            linprogoptions = optimoptions('linprog','Display','none');
            for cdpt_end_i = 1:2*D
                x_end = db_cdpt_end(1:D,cdpt_end_i);
                if any(A_iq*x_end > b_iq,"all") && any(x_end ~= x_n(:,j))
                    d = (x_end > x_n(:,j)) - (x_end < x_n(:,j)); % This is the direction vector from x_n to x_end
                    a = linprog(-1,A_iq*d,b_iq-A_iq*x_n(:,j),[],[],0,1,linprogoptions);
                    db_cdpt_end(1:D,cdpt_end_i) = x_n(:,j) + a*d;
                end
            end
        end
    
        % Append ALL samples that are not sample j, even those within the same batch of new samples
        idx_others    = true(1,X0_len);
        idx_others(j) = false;
        db_cdpt_end = [db_cdpt_end [X0(1:D, idx_others); zeros(CDPT_F_ROWS + CDPT_G_ROWS*g_len + 2, X_n_len-1)]];
    
        % Making sure that all end points are within the (normalized) bounds
        db_cdpt_end = min(1.0, max(0.0, db_cdpt_end));
    
        % Draw candidate points along the directions from x_n to the cdpt_ends
        for cdpt_w = (1:(B-1)) / B
            db_cdpt_iter = [db_cdpt_iter (cdpt_w*db_cdpt_end + (1-cdpt_w)*([X0(1:D,j); zeros(CDPT_F_ROWS + CDPT_G_ROWS*g_len + 2, 1)]))];
        end
    end
else
    db_cdpt_iter    = [];
end

for j = 1:x_n_len
    % Calculate end points from x_n in the D-h_len cardinal directions
    % (ONLY along the axes of the INPUT variables, i.e., not contextual)
    db_cdpt_end = [repmat(x_n(:,j), 1, 2*(D-h_len));
                   zeros(CDPT_F_ROWS + CDPT_G_ROWS*g_len + 2, 2*(D-h_len))];
    n = size(db_cdpt_end, 1);
    db_cdpt_end((h_len+1):(2*n+1):end)     = 0;
    db_cdpt_end((h_len+n+1):(2*n+1):end) = 1;

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
    db_cdpt(1:h_len,:) = repmat(h_n, 1, cdpts);
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