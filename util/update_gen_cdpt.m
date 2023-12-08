%% ==============================
% New candidate points generation
% ===============================

% Build the initial value of the candidate points database chunk to be appended
if isempty(db_cdpt)
    db_cdpt_iter    = [sbl_seq; zeros(CDPT_F_ROWS + CDPT_G_ROWS*g_len + 2, sbl_size)];
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
    % Calculate end points from x_n in the cardinal directions
    db_cdpt_end = [repmat(x_n(:,j), 1, 2*D);
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

new_cdpts = size(db_cdpt_iter, 2);

% Declare new entry set for the cone-generator database
db_cdpt_cn_iter = zeros(2+2*g_len, new_cdpts);

for new_cdpt = 1:new_cdpts

    % Calculate the upper and lower bounds w.r.t. f
    % - *_ulb_lookup stores the locations of the points x_n, and the
    %   heights of the cones from x_n's
    pt_dst = vecnorm(X_n(1:D, :) - repmat(db_cdpt_iter(1:D, new_cdpt), 1, X_n_len));
    f_ulb_lookup = [X_n(XN_FVAL, :); fgam*pt_dst];
    [~, ub_idx] = min(f_ulb_lookup(1, :) + f_ulb_lookup(2, :) + T_n(2, :)); % Calculating individual upper bound values TODO: should add the time-based quantities here
    [~, lb_idx] = max(f_ulb_lookup(1, :) - f_ulb_lookup(2, :) - T_n(2, :)); % Calculating individual lower bound values        
    db_cdpt_iter(CDPT_F_UB_VTX:CDPT_F_LB_HGT, new_cdpt) = [f_ulb_lookup(1:2, ub_idx);
                                                                   f_ulb_lookup(1:2, lb_idx)];
    db_cdpt_cn_iter([CDPTCN_UB CDPTCN_LB], new_cdpt) = [X_n(XN_ID, ub_idx) X_n(XN_ID, lb_idx)]'; % Populating with the index of the cone generator samples

    % Calculate the upper and lower bounds w.r.t. g
    for g_i = 1:g_len
        XN_GVAL_I = XN_GVAL + g_i - 1;
        CDPT_G_I      = CDPT_G_INFO + CDPT_G_ROWS*(g_i-1);
        g_ulb_lookup = [X_n(XN_GVAL_I, :); ggam(g_i)*pt_dst];
        [~, ub_idx] = min(g_ulb_lookup(1, :) + g_ulb_lookup(2, :) + T_n(2+g_i, :)); % Calculating individual upper bound values
        [~, lb_idx] = max(g_ulb_lookup(1, :) - g_ulb_lookup(2, :) - T_n(2+g_i, :)); % Calculating individual lower bound values
        db_cdpt_iter(CDPT_G_I + (CDPT_G_UB_VTX:CDPT_G_LB_HGT), new_cdpt) = [g_ulb_lookup(1:2, ub_idx);
                                                                            g_ulb_lookup(1:2, lb_idx)];
        db_cdpt_cn_iter(2*g_i + [CDPTCN_UB CDPTCN_LB], new_cdpt) = [X_n(XN_ID, ub_idx) X_n(XN_ID, lb_idx)]'; % Populating with the index of the cone generator samples
    end
    db_cdpt_iter(end-1, new_cdpt) = min(vecnorm(X_n(1:D, :) - repmat(db_cdpt_iter(1:D, new_cdpt), 1, X_n_len))) / diam;
end

% Append new candidate points to the big database
db_cdpt    = [db_cdpt db_cdpt_iter];
% Append new data for the cone-generator table to the corresponding big database
db_cdpt_cn = [db_cdpt_cn db_cdpt_cn_iter];

% Calculating the same bounds for the incoming NEW SAMPLE
% the cone lookup should already calculate the time-based uncertainty things

% Generating db_samp entries
if (X_n_len > 1) && timev
    if ~exist('db_samp', 'var')
        db_samp = [];
    end

    if ~exist('db_samp_cn', 'var')
        db_samp_cn = [];
    end

    if isempty(db_samp)
        db_samp_iter    = NaN(4+4*g_len, X_n_len); %
        db_samp_cn_iter = NaN(2+2*g_len, X_n_len); %
        db_samp_ctr     = X_n_len;
    else
        db_samp_iter    = NaN(4+4*g_len, x_n_len); %
        db_samp_cn_iter = NaN(2+2*g_len, x_n_len); %
        db_samp_ctr     = x_n_len;
    end

    for j = 1:db_samp_ctr
        idx_others    = true(1,X_n_len);
        if isempty(db_samp)
            idx_others(j) = false;
            pt_dst       = vecnorm(X_n(1:D, idx_others) - repmat(X_n(1:D,j), 1, X_n_len-1)); %
        else
            idx_others(X_n_len_old+j) = false;
            pt_dst       = vecnorm(X_n(1:D, idx_others) - repmat(x_n(:,j), 1, X_n_len-1));
        end

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

    db_samp    = [db_samp db_samp_iter];
    db_samp_cn = [db_samp_cn db_samp_cn_iter];
end

% elseif X_n_len == 2
%     db_samp = NaN(4+4*g_len, 2);
%     pt_dst = norm(X_n(1:D, 1) - X_n(1:D, 2));
%     db_samp(SAMP_F_UB_VTX:SAMP_F_LB_HGT, :) = [X_n(XN_FVAL, 2) X_n(XN_FVAL, 1);
%                                                        fgam*pt_dst         fgam*pt_dst;
%                                                        X_n(XN_FVAL, 2) X_n(XN_FVAL, 1);
%                                                        fgam*pt_dst         fgam*pt_dst];
%     for g_i = 1:g_len
%         XN_GVAL_I = XN_GVAL + g_i - 1;
%         SAMP_G_I      = SAMP_G_INFO + 4*(g_i-1);
%         db_samp(SAMP_G_I + (SAMP_G_UB_VTX:SAMP_G_LB_HGT), :) = [X_n(XN_GVAL_I, 2) X_n(XN_GVAL_I, 1);
%                                                                 fgam*pt_dst           fgam*pt_dst;
%                                                                 X_n(XN_GVAL_I, 2) X_n(XN_GVAL_I, 1);
%                                                                 fgam*pt_dst           fgam*pt_dst];
% 
%     end
% 
%     % db_samp_cn is a [(2 + 2*g_len) x n] matrix containing the sample index (ID)
%     % of the sample that generated their bounds. Each column represents a sample
%     % (corresponds to each column in X_n). The row enumeration is:
%     % - row 1: the index of the sample which generated the objective upper-bound for this sample
%     % - row 2: the index of the sample which generated the objective lower-bound for this sample
%     % - rows (3,4) and beyond: row pairs similar to rows 1-2, but for each constraint
%     db_samp_cn = repmat([2 1], 2+2*g_len, 1); 
% end

% Calculating the respective uncertainties for the objective and the constraints, storing them in the assigned rows
db_cdpt(CDPT_F_LAMBDA, :) = (db_cdpt(CDPT_F_UB_VTX, :) + db_cdpt(CDPT_F_UB_HGT, :)) ...          % lambda is upper bound ...
                              - (db_cdpt(CDPT_F_LB_VTX, :) - db_cdpt(CDPT_F_LB_HGT, :)) + 2*feps;    % minus the lower bound    
for g_i = 1:g_len
    CDPT_G_I = CDPT_G_INFO + CDPT_G_ROWS*(g_i-1);
    % Uncertainty for constraint
    db_cdpt(CDPT_G_I + CDPT_G_LAMBDA, :) = ...
        (db_cdpt(CDPT_G_I + CDPT_G_UB_VTX, :) + db_cdpt(CDPT_G_I + CDPT_G_UB_HGT, :)) ... 
      - (db_cdpt(CDPT_G_I + CDPT_G_LB_VTX, :) - db_cdpt(CDPT_G_I + CDPT_G_LB_HGT, :)) + 2*geps(g_i);
    % Estimate for constraint (used for predicting if candidate point is feasible or not)
    db_cdpt(CDPT_G_I + CDPT_G_EST, :)    = ...
        sdelta*((db_cdpt(CDPT_G_I + CDPT_G_UB_VTX, :) + db_cdpt(CDPT_G_I + CDPT_G_UB_HGT, :)) ... 
               + (db_cdpt(CDPT_G_I + CDPT_G_LB_VTX, :) - db_cdpt(CDPT_G_I + CDPT_G_LB_HGT, :))) / 2 + ...
        (1 - sdelta)*(db_cdpt(CDPT_G_I + CDPT_G_UB_VTX, :) + db_cdpt(CDPT_G_I + CDPT_G_UB_HGT, :));
end