% ==========================================
% Set Membership Global Optimization (SMGO)
% ==========================================
% minimum usage: smgo(options)
% Necessary arguments to this function (contained in struct 'options'):
%   .objfun  - test function/experiment, accepts a function handle
%   .bounds  - lower and upper bounds. must be an (D x 2) matrix
%   .maxiter - maximum iterations for the optimization run
% For more advanced information and use examples, please refer to the 'Getting Started' guide included with the toolbox.

function out = smgo(options)

defines_const;

%% Acquiring SMGO arguments/options from 'options' struct

get_args;
if err_flag, out = []; return; end

diam = sqrt(D); % diameter of search space (included are decision and context variables), 
                % used in update_gen_cdpt.m and update_db_cdpt.m

%% =============================
% Define initial variable values
% ==============================

% Lipschitz constants
fgam_0    = 1e-6;   ggam_0    = 1e-6*ones(g_len, 1);
fgam      = fgam_0; ggam      = ggam_0;
fgam_prev = fgam_0; ggam_prev = ggam_0; % Used in update_db_cdpt.m

% Noise estimates
feps_0    = 0;      geps_0    = zeros(g_len, 1);
feps      = feps_0; geps      = geps_0;

% Trust region size
tr_exp_0  = 0;
tr_exp    = tr_exp_0;

% Time-varying-related variables
if timev
    % dfeps_dt (dgeps_dt) are quantities that we add to the uncertainties of
    % the objective (constraints) at old samples when the incoming sample 
    % violates the respective existing SM-based bounds.
    dfeps_dt_0 = 1e-6;       dgeps_dt_0 = 1e-6*ones(g_len, 1);
    dfeps_dt   = dfeps_dt_0; dgeps_dt   = dgeps_dt_0;

    % x_del_mdl (x_del_age) hold the samples which are discarded due to
    % model (un)informativeness (too old age)
    x_del_mdl  = [];         x_del_age  = [];
end

% Definition of big databases of samples
% - X_n is a [(D + 2 + S) x N] matrix
% - each column is a valid sample, i.e., not yet discarded
% - for each column
%   - rows 1 to D are the location of the sampling point
%   - row D + 1 is the objective sample
%   - rows D + 2 and beyond are the constraint samples (if there are constraints)
%   - END ROW contains the sample index/ID (currently tied to the iteration in which it was evaluated)
% TODO: decouple the iteration with the sample index (for parallel/asynchronous evaluations)
X_n      = [];    % database of samples (only contains valid samples, i.e., not yet discarded)

% Database of stacked uncertainties (used in the time-varying case)
% - T_n is a (2 + S) x N matrix
% - each column is a valid sample (corresponds to a column in X_n)
% - for each column
%   - row 1 is the sample ID
%   - row 2 is the stacked uncertainty for objective
%   - rows 3 and beyond are stacked uncertainties for the constraints (if present)
T_n      = [];

% Simple database of assigned trust region sizes (used in the contextual case)
% - each column corresponds to a column in X_n
% - contains the EXPONENTS to the trust region size coefficient
%   (0 means default size, higher means smaller trust region, lower means larger trust region)
% NOTE: this database is updated/appended in Algorithm part (3)
R_n      = [];

% Some more databases (for time-varying case)
X_n_all  = []; % Database of ALL sampled points (all valid, and ALSO discarded samples)
X_lf_all = []; % Database of sampled points lifetimes, i.e., the respective ages when discarded

% ============================================
% Definition of output variables and histories
% ============================================
opt_x       = NaN(D, 1);            hist_opt_x  = NaN(D, max_iter);
opt_z       = Inf;                  hist_opt_z  = NaN(1, max_iter);
opt_c       = NaN(g_len, 1);        hist_opt_c  = NaN(g_len, max_iter);

hist_t       = NaN(1, max_iter);
hist_eps     = NaN(1 + g_len, max_iter);
hist_gam     = NaN(1 + g_len, max_iter);
hist_tr      = NaN(1, max_iter); 
hist_deps_dt = NaN(1 + g_len, max_iter);
hist_mode    = NaN(1, max_iter);
mode_prev    = MODE_EXPLOIT; % Used in update_tr_exp.m
% ===========================================

%% Construct the database of segments and candidate points

defines_db;

% Declaration of db_cdpt: from existing points, we calculate the candidate 
%   points (cdpts), and the value of uncertainty. Starting as empty
db_cdpt    = []; % Stores the SM-based upper- and lower-bounds at the candidate points
db_cdpt_cn = []; % Stores the index of the point that generated the (upper-/lower bound) cones
db_samp    = []; % Stores the upper- and lower-bounds information AT the samples (DO NOT DELETE: variable is used in update_rem_samp.m)

%% ============================================================
% Main SMGO loop
% =============================================================
% (1) Evaluation of (black-box) objective and constraint functions
% (2) Calculation of noise bounds and Lipschitz constants, and (only for time-varying case) violation to existing bounds
% (3) Update of the best sample
% (4) Set membership (SM) model building (the bulk part of the algorithm in terms of memory and calc)
% (4.1) (Iterative) update of SM bounds information on existing candidate points
% (4.2) Generation of new candidate points and calculating their SM bounds
% (5) Exploitation routine (simply select a candidate point, using the updated SM bounds)
% (6) Exploration routine (similarly just a selection process of candidate point)

for iter = 1:max_iter
    if iter == 1      
        % If initial samples set is already available, put into database already
        % X0 samples SHOULD already include the corresponding contexts (if applicable)
        if exist('X0', 'var')
            X0_nrm  = real2normd(X0(1:D,:), bnds);
            X0_len  = size(X0_nrm, 2);
            X_n     = [X0_nrm; X0((D+1):end,:); 1:X0_len; zeros(1,X0_len)];
            X_n_all = [X0_nrm; X0((D+1):end,:)];
            T_n     = [1:X0_len; zeros(1 + g_len, X0_len)];
            R_n     = zeros(1,X0_len);

            xn_id   = X0_len + 1;
        else
            X0_len  = 0; 
        end

        % Assign the initial sampling point
        x_n = x0;  
        
        % If with context, evaluate the context NOW, and 
        % associate with the initial sampling point
        if h_len % || exact_calc
            % Get all context values
            h_n = NaN(h_len,1);
            if isfield(options, 'cxtfun')
                % for h_i = 1:h_len
                    % Context variables are NOT functions of the input x_n,
                    % these are simply measurements of the environment
                    h_n = h{1}(); % Definition of h in get_args.m as list of function handles
                % end
                h_n = real2normd(h_n,bnds(1:h_len,:));
            end
            if ~isempty(x_n)
                x_n(1:h_len,:) = repmat(h_n,1,size(x_n,2));
            end
        end
    end
       
    %% =================
    % Algorithm part (1)
    % ==================

    % BEFORE evaluating and appending the new sample/s (with zero time-based uncertainty) to T_n, 
    % I shall stack the previously-calculated uncertainty violations to the existing entries
    if timev
        if ~isempty(T_n)
            T_n(2:end, :) = T_n(2:end, :) + repmat([dfeps_dt; dgeps_dt], 1, size(T_n,2));
        end
    end
    
    % Tracking the number of entries in X_n before the introduction of new samples
    X_n_len_old = size(X_n, 2);
    % Had the measure the number of columns because it can be that x_n represents multiple sampling points
    x_n_len     = size(x_n, 2);
    
    % Iterate over each sampling point in x_n, and evaluate the objective and constraints
    for j = 1:x_n_len
        zc = NaN(1 + g_len, 1);

        % If TCP interface has been set up
        if tcp_en
            % Read the data received by the SMGO TCP server
            % Data should be a string composed of space-separated numerical values
            % (American convention, period as decimal separator)
            % - objective value
            % - constraint value/s (if present)
            % - context variable measurement/s (if present)
            snd_dat = num2str(normd2real(x_n,bnds)');
            fprintf("Sent %s\n", snd_dat);
            write(f_tcp,snd_dat,"string");
            while ~f_tcp.NumBytesAvailable
                pause(0.01);
            end
            rcv_dat = read(f_tcp,f_tcp.NumBytesAvailable,"string");
            fprintf("Received %s\n",rcv_dat);
            zc = str2num(rcv_dat);  
        else
            if ~isfield(options, 'confun') % && ~isfield(options, 'cxtfun')
                % If the objective and constraints values are already
                % contained in the 'f' function as a vector
                if ~isempty(f)
                    zc = f(normd2real(x_n(:,j), bnds)); 
                else
                    zc = [];
                end
            else
                % If the constraint functions are enumerated separately
                % in the field options.confun. NOTE: if h_len > 0 (contextual
                % case), x_n ALREADY contains the context!!! See Algorithm (4)b
                zc(1) = f(normd2real(x_n(:,j), bnds));

                % Run all constraint functions
                if isfield(options, 'confun')
                    for g_i = 1:g_len
                        zc(1 + g_i) = g{g_i}(normd2real(x_n(:,j), bnds));
                    end  
                end
            end
        end
        
        if ~isempty(zc)
            if size(zc,2) ~= 1
                zc = zc';
            end

            z_n = zc(1);
            if g_len > 0
                c_n = zc(2:(1+g_len));
            else
                c_n = [];
            end

            % Stacking the new sample to the big collection
            % NOTE AGAIN: x_n ALREADY contains the context (please see Algorithm (4)b to (6))
            X_n     = [X_n [x_n(:,j); z_n; c_n; xn_id; iter]];   % Side-stacking the column vector to collection of samples (fat matrix)
            X_n_all = [X_n_all [x_n(:,j); z_n; c_n]];
            T_n     = [T_n [xn_id; zeros(1 + g_len, 1)]];        % Similarly stacking the time-related column vector to the collection 
                                                                 % TODO: 'iter' should be the timestamp instead 
                                                                 % TODO: why do I have the zeros** vector here? ANSWER: see definition for T_n (smgo.m, lines 63-69)
            R_n     = [R_n NaN];                                 % Stacking the TR radius vector with a dummy value (will be updated with a real one in Algorithm part (3))
            xn_id   = xn_id + 1;
        end
    end
    
    X_n_len = size(X_n, 2);

    % =========================
    % End of algorithm part (1)
    % =========================

    %% =================
    % Algorithm part (2)
    % ==================
    calc_time = tic;

    if timev && (iter == 1)
        calc_eps; % Estimate noise bounds
        calc_gam; % Estimate Lipschitz constants
    end

    if iter > 1
        if ~timev
            calc_eps; % Estimate noise bounds
            calc_gam; % Estimate Lipschitz constants
        end

        % Time-varying case: calculation of degree of violation of existing bounds
        % TODO: STILL NOT accommodating context variables
        if timev
            dfeps_dt = 0.0;
            for j = 1:x_n_len
                pt_dst = vecnorm(X_n(1:D, 1:X_n_len_old) - repmat(x_n(:,j), 1, X_n_len_old));
                f_up_bnd = min(X_n(XN_FVAL, 1:X_n_len_old) + fgam_prev*pt_dst); % the existing (old) upper bound at the location of new sample
                f_lo_bnd = max(X_n(XN_FVAL, 1:X_n_len_old) - fgam_prev*pt_dst); % the existing (old) lower bound
                
                dfeps_up = max((X_n(XN_FVAL, X_n_len_old+j) + feps) - f_up_bnd, 0); % if (X_n(XN_FVAL, end) + feps) VIOLATES old upper bound            
                dfeps_lo = max(f_lo_bnd - (X_n(XN_FVAL, X_n_len_old+j) - feps), 0); % if (X_n(XN_FVAL, end) - feps) VIOLATES old lower bound
                
                dfeps_dt = max([dfeps_dt dfeps_up dfeps_lo]);
            end
        end    
    end    
    hist_eps(:, iter) = [feps; geps];
    hist_gam(:, iter) = [fgam; ggam];

    if timev
        hist_deps_dt(:, iter) = [dfeps_dt; dgeps_dt];
    end
    % =========================
    % End of algorithm part (2)
    % =========================

    %% =================
    % Algorithm part (3)
    % ==================

    % Identifying the best feasible sample
    % There are different approaches, the one in the below
    % clause is for the case WITHOUT contextual variables
    if ~h_len
        % Find the minimum value of z among valid samples
        X_n_vld   = X_n(:,all(X_n(XN_GVAL:end-2,:)<=0, 1));
        [opt_z_new, opt_z_idx] = min(X_n_vld(XN_FVAL,:));
        opt_x_new = X_n_vld(1:D,opt_z_idx);
        opt_c_new = X_n_vld(XN_GVAL:end-2,opt_z_idx);
    
        % Updating the current best sample if valid best sample exists
        if (~isempty(opt_z_new) && ~isnan(opt_z_new)) && ...
           ((timev && (opt_z_new ~= opt_z_old)) || ...
            (~timev && (opt_z_new < opt_z_old)))
            opt_z = opt_z_new;
            opt_x = opt_x_new;
            opt_c = opt_c_new;
            opt_z_upd = true; % This variable IS used: in update_tr_exp.m
        else
            opt_z_upd = false;
        end

    % The below clause describes the case WITH context
    else
        % Filter samples which are feasible
        X_n_vld     = X_n(:,all(X_n(XN_GVAL:end-2,1:end-1)>0, 1));
        X_n_vld_len = size(X_n_vld,2);

        % Filter samples which are close
        if h_len == 1 % If context is 1-dimensional
            cxt_dist = abs(X_n_vld(1:h_len,:) - repmat(h_n, 1, X_n_vld_len));
        else          % If context is more-dimensional
            cxt_dist = vecnorm(X_n_vld(1:h_len,:) - repmat(h_n, 1, X_n_vld_len));
        end
        cxt_near = cxt_dist < cxt_rad;

        X_n_near = [];
        if ~isempty(X_n_vld)
            X_n_near = X_n_vld(:, cxt_near);
        end

        % IF X_n_near EXISTS
        if ~isempty(X_n_near)
            [~, opt_z_idx] = min(X_n_near(XN_FVAL,:));
            opt_x = X_n_near(1:D, opt_z_idx); % This opt_x STILL contains its (only nearby) sampled context
                                              % NOT the presently-measured context
            
            % Project data from (opt_x, opt_z) to the just-measured context
            opt_x(1:h_len) = h_n; % Project opt_x to new context

            % TODO: calculate the exact central approximation at the projected opt_x
            pt_dst       = vecnorm(X_n(1:D, :) - repmat(opt_x, 1, X_n_len));
            f_ulb_lookup = [X_n(XN_FVAL, :); fgam*pt_dst];
            opt_x_ub = min(f_ulb_lookup(1, :) + f_ulb_lookup(2, :) + feps + T_n(2, :)); % Calculating individual upper bound values
            opt_x_lb = max(f_ulb_lookup(1, :) - f_ulb_lookup(2, :) - feps - T_n(2, :)); % Calculating individual lower bound values 

            opt_z = (opt_x_ub + opt_x_lb)/2;
            if opt_z < opt_z_old
                opt_z_upd = true;
            else
                opt_z_upd = false;
            end

            % TODO: opt_c MIGHT be unfeasible because of the projections 
            % of the constraint functions
            opt_c = NaN(g_len,1);
            for g_i = 1:g_len
                XN_GVAL_I = XN_GVAL + g_i - 1;
                g_ulb_lookup = [X_n(XN_GVAL_I, :); ggam(g_i)*pt_dst];
                c_ub = min(g_ulb_lookup(1, :) + g_ulb_lookup(2, :) + geps(g_i) + T_n(2+g_i, :));
                c_lb = max(g_ulb_lookup(1, :) - g_ulb_lookup(2, :) - geps(g_i) - T_n(2+g_i, :));
                opt_c(g_i) = (c_ub + c_lb)/2;
            end
            
            % Update/append the TR radius database
            % Trust region dynamics (contextual optimization) happens here
            % - get the sample index of the best (context-related) sample
            % - if generated through exploitation, compare with the best sample
            %   - if better at least by the expected improvement, best sample TR entry decreases
            %   - if NOT better, best sample TR entry increases
            %   - else (not shown anymore in code), best sample TR entry does NOT change
            % - COPY TR entry of new sample from that of best sample
            R_n_idx = X_n_near(XN_ID,opt_z_idx);
            if iter == 1
                R_n(end) = 0;
            else
                if mode_prev == MODE_EXPLOIT
                    if z_n < exploit_thr
                        R_n(R_n_idx) = R_n(R_n_idx) - 1;
                    elseif z_n > opt_z
                        R_n(R_n_idx) = R_n(R_n_idx) + 1;
                    end
                end
                R_n(end) = R_n(R_n_idx);
            end

        % IF X_n_near DOES NOT EXIST
        % then go directly to explore (Algorithm part (6))!
        else
            opt_x = NaN(D, 1);
            opt_z = Inf;
            opt_c = NaN(g_len,1);
            opt_z_upd = false;
            R_n(end) = 0;
        end
    end
    % =========================
    % End of algorithm part (3)
    % =========================

    %% =====================================
    % Algorithm part (4)a: iterative updates
    % ======================================
    % ONLY runs if there are no context variables

    if ~h_len
        % ================================================
        % Iteratively update the candidate points database db_cdpt
        % ================================================
        % - Saves the vertex value of the upper (lower) bound
        % - Also saves the hypercone height value
        % - Only the height value is scaled
        % - If the new incoming upper (lower) bound is tighter than existing
        %   - Replace the vertex value
        %   - Replace the height value
        if iter > 1
            update_db_cdpt;
            % The condition (iter > 2) is here because at iter == 2,
            % I pre-populate the initial contents of db_samp (see update_gen_cdpt.m)
            if timev && (iter > 2)
                update_db_samp;
            end
        end 
        
        % =====================================
        % Introduce additional cdpts to db_cdpt
        % =====================================
        update_gen_cdpt;
        % TODO: update_gen_cdpt also in contextual optimization
        %       - when initial samples are supplied BUT no initial point
    
        % ===========================================================================
        % Remove non-informative samples (only applicable for time-varying case)
        % ===========================================================================
        % TODO: should I remove samples BEFORE generating the candidate points? (see previous code chunk)
        if timev
            update_rem_samp;
            % [2023-01-03] TODO: no need to double check. Scans current contents of the databases
        end
        
        X_n_len = size(X_n, 2); % Update the number of elements in samples list
    end
    % ==========================
    % End of algorithm part (4)a
    % ==========================

    %% ========================================================================
    % Algorithm part (4)b: Measure context and re-calculate all SM-based bounds
    % =========================================================================
    % Runs IF there are contextual variables
    skip2end = false;
    if iter == options.maxiter
        skip2end = true;
    end
    if h_len && ~skip2end% || exact_calc
        % Get all context values
        h_n      = NaN(h_len,1);
        if isfield(options, 'cxtfun')
            % for h_i = 1:h_len
                % Context variables are NOT functions of the input x_n,
                % these are simply measurements of the environment
                h_n = h{1}(); % Definition of h in get_args.m as list of function handles
            % end
            h_n = real2normd(h_n,bnds(1:h_len,:));
        end

        % Regenerate the candidate points, projecting them to the new context
        % and calculating their corresponding bounds from scratch
        update_regen_cdpt;

        % Filter samples which are feasible
        X_n_vld     = X_n(:,all(X_n(XN_GVAL:end-2,1:end-1)>0, 1));
        X_n_vld_len = size(X_n_vld,2);

        % Filter samples which are close
        if h_len == 1 % If context is 1-dimensional
            cxt_dist = abs(X_n_vld(1:h_len,:) - repmat(h_n, 1, X_n_vld_len));
        else          % If context is more-dimensional
            cxt_dist = vecnorm(X_n_vld(1:h_len,:) - repmat(h_n, 1, X_n_vld_len));
        end
        cxt_near = cxt_dist < cxt_rad;

        X_n_near = [];
        if ~isempty(X_n_vld)
            X_n_near = X_n_vld(:, cxt_near);
        end

        % IF X_n_near EXISTS
        if ~isempty(X_n_near)
            [~, opt_z_idx] = min(X_n_near(XN_FVAL,:));
            opt_x = X_n_near(1:D, opt_z_idx); % This opt_x STILL contains its (only nearby) sampled context
                                              % NOT the presently-measured context
            
            % Project data from (opt_x, opt_z) to the just-measured context
            opt_x(1:h_len) = h_n; % Project opt_x to new context

            % TODO: calculate the exact central approximation at the projected opt_x
            pt_dst       = vecnorm(X_n(1:D, :) - repmat(opt_x, 1, X_n_len));
            f_ulb_lookup = [X_n(XN_FVAL, :); fgam*pt_dst];
            opt_x_ub = min(f_ulb_lookup(1, :) + f_ulb_lookup(2, :) + feps + T_n(2, :)); % Calculating individual upper bound values
            opt_x_lb = max(f_ulb_lookup(1, :) - f_ulb_lookup(2, :) - feps - T_n(2, :)); % Calculating individual lower bound values 

            opt_z = (opt_x_ub + opt_x_lb)/2;

            % TODO: opt_c MIGHT be unfeasible because of the projections 
            % of the constraint functions
            opt_c = NaN(g_len,1);
            for g_i = 1:g_len
                XN_GVAL_I = XN_GVAL + g_i - 1;
                g_ulb_lookup = [X_n(XN_GVAL_I, :); ggam(g_i)*pt_dst];
                c_ub = min(g_ulb_lookup(1, :) + g_ulb_lookup(2, :) + geps(g_i) + T_n(2+g_i, :));
                c_lb = max(g_ulb_lookup(1, :) - g_ulb_lookup(2, :) - geps(g_i) - T_n(2+g_i, :));
                opt_c(g_i) = (c_ub + c_lb)/2;
            end
            tr_exp = R_n(opt_z_idx);

        % IF X_n_near DOES NOT EXIST
        % then go directly to explore (Algorithm part (6))!
        else
            opt_x = NaN(D, 1);
            opt_z = Inf;
            opt_c = NaN(g_len,1);
            opt_z_upd = false;
            R_n(end) = 0;
        end
    end

    %% =========================================================
    % Algorithm part (4)c: Generate/update trust region hyperbox
    % ==========================================================

    % IF opt_x exists, then there IS a trust region (and consequently, allow exploitation)
    % - IF difference of current context from previous one is not much,
    %   THEN do a trust region update just like in classical (contextless) SMGO
    % - IF current context is super different, THEN trust region is reset to the default size
    % - EITHER WAY, initialize the Sobol-generated candidate points
    % IF opt_x does NOT exist, then there is NO trust region, and proceed directly to exploration
    if opt_z < Inf && ~skip2end
        update_tr_exp;
        hist_tr(iter) = tr_exp;
        
        if sbl_size > 0
            % Create candidate points inside tr_bounds, decided by location of Sobol set points
            % Matrix sbl_seq is [D x sbl_seq]
            sbl_cdpt = [[sbl_seq*tr_size*(tr_coeff ^ tr_exp); zeros(h_len, sbl_size)] + tr_bnds(:, 1)*ones(1, sbl_size);  % Matrix tr_bnds is defined in update_tr_exp.m
                        zeros(2 + 2*g_len, sbl_size)];
            sbl_cdpt(1:h_len,:) = repmat(h_n, 1, sbl_size);

            % For each Sobol candidate point, calculate SM-bounds from scratch
            for sbl_i = 1:sbl_size
                sbl_dst = vecnorm(X_n(1:D,:) - repmat(sbl_cdpt(1:D, sbl_i), 1, X_n_len));
                sbl_cdpt(SBL_FUB, sbl_i) = min(X_n(XN_FVAL,:) + feps + fgam*sbl_dst + T_n(2, :));  
                sbl_cdpt(SBL_FLB, sbl_i) = max(X_n(XN_FVAL,:) - feps - fgam*sbl_dst - T_n(2, :));
                for g_i = 1:g_len
                    SBL_G_I = SBL_G_INFO + 2*(g_i-1);
                    XN_GVAL_I = XN_GVAL + g_i - 1;
                    sbl_cdpt(SBL_G_I + SBL_GUB, sbl_i) = min(X_n(XN_GVAL_I,:) + geps(g_i) + ggam(g_i)*sbl_dst + T_n(2+g_i, :));  
                    sbl_cdpt(SBL_G_I + SBL_GLB, sbl_i) = max(X_n(XN_GVAL_I,:) - geps(g_i) - ggam(g_i)*sbl_dst - T_n(2+g_i, :));
                end
            end
        end
    end

    %% ================================================
    % Algorithm part (5): sampling mode 1, exploitation
    % =================================================
    
    % The exploitation should remain the same, and will not consider the time-related uncertainty. 
    % It should remain only dependent on the uncertainty w.r.t. existing (valid) samples.
    mode1_ok = false;
    if opt_z < Inf && ~skip2end
        % Processing and choosing exploitation point from candidate points database
        cdpts        = size(db_cdpt, 2);
        cdpt_vld_idx = true(1, cdpts);
        cdpt_vld_idx = cdpt_vld_idx & all(db_cdpt(1:D,:) >= 0) & all(db_cdpt(1:D,:) <= 1);
    
        % Filtering due to constraints satisfaction
        for g_i = 1:g_len
            CDPT_G_I     = CDPT_G_INFO + CDPT_G_ROWS*(g_i-1);
            cdpt_vld_idx = cdpt_vld_idx & (db_cdpt(CDPT_G_I + CDPT_G_EST, :) <= 0);
        end    
    
        % Filtering due to trust region
        cdpt_vld_idx = cdpt_vld_idx & prod(db_cdpt(1:D, :) >= tr_bnds(:, 1)); % The product ('prod') operation is performed PER COLUMN.
        cdpt_vld_idx = cdpt_vld_idx & prod(db_cdpt(1:D, :) <= tr_bnds(:, 2)); % The resulting (row) vector is AND-gated
        
        % Filtering due to linear inequalities (if present)
        if isfield(options,'ineq')
            cdpt_vld_idx = cdpt_vld_idx & all(A_iq*db_cdpt(1:D, :) <= repmat(b_iq,1,cdpts));
        end

        % TODO: avoid making copies of filtered out matrices
        if sum(cdpt_vld_idx) > 0
            if ~h_len
                w_min = ((db_cdpt(CDPT_F_UB_VTX, :) + db_cdpt(CDPT_F_UB_HGT, :)) + ...
                         (db_cdpt(CDPT_F_LB_VTX, :) - db_cdpt(CDPT_F_LB_HGT, :))) / 2;
                w_vld = -w_min + sbeta*db_cdpt(CDPT_F_LAMBDA, :);
                [~, cdpt_idx] = max(w_vld + (~cdpt_vld_idx*-INFTY));
            else
                w_min = (db_cdpt(CDPT_F_UB, :) + db_cdpt(CDPT_F_LB, :)) / 2;
                w_vld = -w_min + sbeta*db_cdpt(CDPT_F_LAMBDA, :);
                [~, cdpt_idx] = max(w_vld + (~cdpt_vld_idx*-INFTY));
            end
    
            % The 'while' snippet below is a last-resort way to avoid sampling
            % at a previously-sampled location, but practically should not run at all!
            cdpt_tmp = db_cdpt(1:D, cdpt_idx);
            % If chosen point x_n is too close to an existing sample, choose another
            while ~isempty(cdpt_tmp) && min(vecnorm(X_n(1:D, :) - repmat(cdpt_tmp, 1, X_n_len))) < 1e-9
    %             fprintf('Iteration %d: Choosing another exploitation point. Distance is %e\n', iter, min(vecnorm(X_n(1:D, :) - repmat(cdpt_tmp, 1, X_n_len))));
                w_vld(:, cdpt_idx)      = [];
                cdpt_vld_idx(cdpt_idx)  = [];
                db_cdpt(:, cdpt_idx)    = [];
                if ~h_len
                    db_cdpt_cn(:, cdpt_idx) = [];
                end
                [~, cdpt_idx]         = max(w_vld + (~cdpt_vld_idx*-INFTY));
                cdpt_tmp              = db_cdpt(1:D, cdpt_idx);
            end
        else
            cdpt_tmp = [];
        end
        
        if sbl_size > 0
            % Processing and choosing exploitation point from Sobol-generated candidate points
            sbl_vld_idx   = true(1, sbl_size);
            sbl_vld_idx   = sbl_vld_idx & all(sbl_cdpt(1:D,:) >= 0) & all(sbl_cdpt(1:D,:) <= 1);
            % Filtering due to constraints satisfaction
            for g_i = 1:g_len
                SBL_G_I = SBL_G_INFO + 2*(g_i-1);
                sbl_vld_idx = sbl_vld_idx & (sdelta*(sbl_cdpt(SBL_G_I + SBL_GUB, :) + sbl_cdpt(SBL_G_I + SBL_GLB, :)) / 2 + ...
                                              (1 - sdelta)*sbl_cdpt(SBL_G_I + SBL_GUB, :) <= 0);
            end

            % Filtering due to linear inequalities (if present)
            if isfield(options,'ineq')
                sbl_vld_idx = sbl_vld_idx & all(A_iq*sbl_cdpt(1:D, :) <= repmat(b_iq,1,sbl_size));
            end
            sbl_vld = sbl_cdpt(:, sbl_vld_idx);
            if ~isempty(sbl_vld)
                [~, cdpt_idx] = max(-(sbl_vld(SBL_FUB, :) + sbl_vld(SBL_FLB, :)) / 2 + sbeta*(sbl_vld(SBL_FUB, :) - sbl_vld(SBL_FLB, :)));
        
                % This 'while' snippet below is a last-resort way to avoid sampling 
                % at a previously-sampled location but practically should not run at all!
                sbl_tmp = sbl_vld(1:D, cdpt_idx);
        
                % If chosen Sobol point x_n is too close to an existing sample, choose another
                while ~isempty(sbl_tmp) && min(vecnorm(X_n(1:D, :) - repmat(sbl_tmp, 1, X_n_len))) < 1e-9         
                    sbl_vld(:, cdpt_idx) = [];
                    [~, cdpt_idx]        = min(sbl_vld(end, :));        
                    sbl_tmp              = sbl_vld(1:D, cdpt_idx);            
                end
            else
                sbl_tmp = [];
            end
        else
            sbl_tmp = [];
        end
        
        % Just selecting if the candidate point from db_cdpt or sbl_cdpt should be chosen
        if ~isempty(cdpt_tmp) && ~isempty(sbl_tmp)
            if max((X_n(XN_FVAL,:) - feps) - fgam*vecnorm(X_n(1:D,:) - repmat(cdpt_tmp, 1, X_n_len))) < ...
               max((X_n(XN_FVAL,:) - feps) - fgam*vecnorm(X_n(1:D,:) - repmat(sbl_tmp, 1, X_n_len)))
                x_n_tmp = cdpt_tmp;
            else
                x_n_tmp = sbl_tmp;
            end
        elseif ~isempty(cdpt_tmp)
            x_n_tmp = cdpt_tmp;        
        elseif ~isempty(sbl_tmp)
            x_n_tmp = sbl_tmp;
        else
            x_n_tmp = [];
        end
        
        exploit_thr = opt_z - feps - salpha*fgam;
        if ~isempty(x_n_tmp)
            if opt_z < Inf && max((X_n(XN_FVAL,:) - feps) - fgam*vecnorm(X_n(1:D,:) - repmat(x_n_tmp, 1, X_n_len))) < exploit_thr
                x_n = x_n_tmp;
                hist_mode(iter) = MODE_EXPLOIT;
                mode_prev       = MODE_EXPLOIT;
                mode1_ok        = true;
            end
        end
    end
    % =========================
    % End of algorithm part (5)
    % =========================
    
    %% ===============================================
    % Algorithm part (6): sampling mode 2, exploration
    % ================================================

    cdpts        = size(db_cdpt, 2);
    cdpt_vld_idx = true(1, cdpts);
    w_bst        = ones(1, cdpts);
    w_unc        = zeros(1, cdpts); 
    
    % =================================
    % Calculating merit function values
    % =================================
    for g_i = 1:g_len
        CDPT_G_I     = CDPT_G_INFO + CDPT_G_ROWS*(g_i-1);
        g_i_vld      = (db_cdpt(CDPT_G_I + CDPT_G_EST, :) <= 0);
        w_bst        = w_bst .* (g_i_vld+1);
        w_unc        = w_unc + (db_cdpt(CDPT_G_I + CDPT_G_LAMBDA, :))/ggam(g_i);
        cdpt_vld_idx = cdpt_vld_idx & g_i_vld;
    end

    w_vld = (db_cdpt(CDPT_F_LAMBDA, :)) .* cdpt_vld_idx / fgam;
    mode2_mrt = ((1-sdelta)*w_vld + sdelta*w_unc.*w_bst/(2^g_len)).*db_cdpt(end-1, :) + phi.*db_cdpt(end, :);
    
    
    % ==========================================================
    % Clean up some cdpts from db_cdpt which are most likely 
    % not selected up until the next filtercdpt_horiz iterations
    % ==========================================================
    % TODO: Should edit for contextual optimization
    %       - the entry should NOT be deleted
    %       - instead, the age should just be reset
    if filtercdpt && (iter >= filtercdpt_horiz) && ~skip2end
        cdpt_clnup_idx = ((max(mode2_mrt) - mode2_mrt)/phi >= filtercdpt_horiz);
        cdpt_clnup_idx(1:sbl_size)    = 0;
        db_cdpt(:, cdpt_clnup_idx)    = [];
        if ~h_len
            db_cdpt_cn(:, cdpt_clnup_idx) = [];
        end
        w_vld(:, cdpt_clnup_idx)      = [];
        w_unc(:, cdpt_clnup_idx)      = [];
        w_bst(:, cdpt_clnup_idx)      = [];
        mode2_mrt = ((1-sdelta)*w_vld + sdelta*w_unc.*w_bst/(2^g_len)).*db_cdpt(end-1, :) + phi.*db_cdpt(end, :);
    end

    % ================================
    % The actual exploration routine
    % ================================
    if ~mode1_ok && ~skip2end
        [~, cdpt_idx] = max(mode2_mrt);
        x_n = db_cdpt(1:D, cdpt_idx);

        % This 'while' snippet below is a last-resort way to avoid sampling 
        % at a previously-sampled location but practically should not run at all!
        while ~isempty(x_n) && min(vecnorm(X_n(1:D, :) - repmat(x_n, 1, X_n_len))) < 1e-9  
            % If chosen point x_n is too close to an existing sample, 
            % delete that candidate point from the database, and choose another
%             fprintf('Iteration %d: Choosing another exploration point. Distance is %e\n', iter, min(vecnorm(X_n(1:D, :) - repmat(x_n, 1, X_n_len))));
            db_cdpt(:, cdpt_idx)    = [];
            if ~h_len
                db_cdpt_cn(:, cdpt_idx) = [];
            end
            w_vld(:, cdpt_idx)      = [];
            w_unc(:, cdpt_idx)      = [];
            w_bst(:, cdpt_idx)      = [];

            mode2_mrt = ((1-sdelta)*w_vld + sdelta*w_unc.*w_bst/(2^g_len)).*db_cdpt(end-1, :) + phi.*db_cdpt(end, :);
            
            [~, cdpt_idx] = max(mode2_mrt);
            x_n = db_cdpt(1:D, cdpt_idx);
        end
        hist_mode(iter) = MODE_EXPLORE;
        mode_prev       = MODE_EXPLORE;
    end
    % =========================
    % End of algorithm part (6)
    % =========================

    % If no context, we can remove the sampling
    % point from the candidate points database
    cdpt_idx = all(db_cdpt(1:D, :) == x_n);
    if ~h_len
        db_cdpt(:, cdpt_idx)    = [];
        db_cdpt_cn(:, cdpt_idx) = [];
        w_vld(:, cdpt_idx)      = [];
        w_unc(:, cdpt_idx)      = [];
        w_bst(:, cdpt_idx)      = [];

    % If instead we have contexts, we DO NOT
    % remove the sampling point from the database,
    % instead we just reset the age to zero
    else
        db_cdpt(end, cdpt_idx)  = 0;
    end

    % Clip x_n inside the (normalized) search space
    x_n = min(1.0, max(0.0, x_n));

    if tcp_en
        write(f_tcp, x_n, "double");
    end

    hist_t(iter) = toc(calc_time);
    if ~exist('opt_x', 'var')
        opt_x = NaN(D, 1);
    end    
    hist_opt_x(:, iter) = opt_x;
    hist_opt_z(iter)    = opt_z;
    hist_opt_c(:, iter) = opt_c;

    feps_prev = feps;
    geps_prev = geps;
    opt_z_old = opt_z;

    % Running per-iteration feedback function
    if ~isempty(iterfdbackfun)
        iterfdbackfun(X_n);
    end
    
    if iter == 1
        fprintf("============\nRunning SMGO\n============\n");
    end
    % fprintf("Iteration: %d/%d, ",iter, max_iter);
    % if ~h_len % Only display a 'best value' when not in contextual optimization. `opt_z` does not make sense in contextual case.
    %     if opt_z < Inf
    %         fprintf("Best z: %.3f\n", opt_z);
    %     else
    %         fprintf("No feasible point found yet.");
    %     end
    % end
    % fprintf("\n");
end
fprintf("Optimization complete. Results are stored inside the output struct.\n")

X_lf_all = [X_lf_all [X_n(1:D,:); iter - X_n(XN_TSTMP,:)]];

hist_x   = X_n(1:D,:);    
hist_z   = X_n(XN_FVAL,:);
hist_c   = X_n(XN_GVAL:end-2,:);
x_all    = X_n_all;
x_lf_all = X_lf_all;

%% ==================
% Build output struct
% ===================
out.opt_x = normd2real(opt_x, bnds);        out.hist_x = normd2real(hist_x, bnds);       out.hist_opt_x = normd2real(hist_opt_x, bnds);
out.opt_z = opt_z;                          out.hist_z = hist_z;                         out.hist_opt_z = hist_opt_z;
out.nxt_x = normd2real(x_n, bnds);          out.hist_x_normd = hist_x;

if g_len > 0
    out.opt_c = opt_c;
    out.hist_c = hist_c;       
    out.hist_opt_c = hist_opt_c;
end
                          
out.hist_t    = hist_t;
out.hist_mode = hist_mode;
out.hist_eps  = hist_eps;
out.hist_gam  = hist_gam;
out.hist_tr   = hist_tr;

if timev
    out.hist_x_all   = x_all;
    out.x_life_t     = x_lf_all;
    out.x_del_mdl    = x_del_mdl;
    out.x_del_age    = x_del_age;
    out.hist_deps_dt = hist_deps_dt;
end

% Running output feedback function
if ~isempty(outfdbackfun)
    outfdbackfun(out);
end

%     if D == 2
%         % for plotting
%         pso_options = optimoptions('particleswarm','Display','none');
%         [xmin_pso, zmin_pso] = particleswarm(@(x)f(x, max(iter-timev_init,0)),2,[-10 -10],[10 10],pso_options);
%     end


%     if D == 2
%         % ===========
%         % plotting
%         % ===========
%     
%         plot3(X_n(1,:)*10-5,X_n(2,:)*10-5,X_n(3,:),'.','MarkerSize',15); hold on;
%         plot(opt_x(1,:)*10-5,opt_x(2,:)*10-5,'x','LineWidth',2,'MarkerSize',12);
%         view(0,90);
%     
%         x_ax = linspace(-5,5,100);
%         y_ax = linspace(-5,5,100);
%         grid_ax = meshgrid(x_ax,y_ax);
%         z_v = zeros(size(grid_ax));
%         for i = 1:length(x_ax)
%             for j = 1:length(y_ax)
%                 z_v(j,i) = f([x_ax(i) y_ax(j)], max(iter-timev_init,0));
%             end
%         end
%     
%         contour(x_ax,y_ax,z_v, 25, 'ShowText', 'off');
%         if exist('del_pts_plt_mdl','var')
%             plot( del_pts_plt_mdl(1,:)*10-5, del_pts_plt_mdl(2,:)*10-5, 'ro');
%         end
%     
%         if exist('del_pts_plt_age','var')
%             plot( del_pts_plt_age(1,:)*10-5, del_pts_plt_age(2,:)*10-5, 'r*');
%         end
% 
% 
%         title(sprintf('Evol = %.3f%%, Min = %.3f, Value @ opt-x = %.3f, Gap = %.3f',min(1, max(iter-timev_init,0)/3000)*100, zmin_pso, z_opt_x, z_opt_x_gap));
%         drawnow;
%         hold off;
%     end
