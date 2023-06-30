% ===========================================
% Extract the absolutely required parameters
% ===========================================
max_iter = options.maxiter;
bnds     = options.bounds;
D        = size(bnds, 1); % D includes the input variables and context variables
% ===========================================

err_flag = false;

% ===========================================
% Objective function parameters
% ===========================================
tcp_en = false;
if isstruct(options.objfun)
    if isfield(options.objfun, 'addr') && isfield(options.objfun, 'port')
        tcp_settings = options.objfun;
        fprintf("[STATUS] Opening SMGO TCP server at %s:%d...\n",tcp_settings.addr,tcp_settings.port);
        pause(0.01);
        f_tcp = tcpserver(tcp_settings.addr, tcp_settings.port);
        % If error occurs in the tcpserver call, let the tcpserver just end the SMGO runtime and return the error.
        fprintf("[STATUS] Waiting for client to connect...\n");
        pause(0.01);
        while ~f_tcp.Connected
            pause(0.01);
        end
        tcp_en = true;
        fprintf("[STATUS] Client connected.\n");
        pause(0.01);
    else
        fprintf("[ERROR] Invalid settings for TCP interface. Please supply options.objfun with fields addr and port.\n");
        err_flag = true;
    end
elseif isa(options.objfun, 'function_handle')
    f = options.objfun;
else
    fprintf("[ERROR] Wrong objfun argument. Please supply with a function handle or a struct with fields addr and port.\n")
    err_flag = true;
end
% ===========================================

% ===========================================
% Black-box constraint function/s parameters
% ===========================================
if isfield(options, 'confun')
    g = options.confun;
    g_len = length(g);
    if isa(g, 'function_handle')
        g = {g};
    end
elseif isfield(options, 'conlen')
    g_len = options.conlen;
else
    g_len = 0;
    g = {};
end
% ===========================================

% ===========================================
% Linear inequality constraints
% ===========================================
if isfield(options, 'ineq')
    A_iq = options.ineq.A * diag(bnds(:,2)-bnds(:,1));
    b_iq = options.ineq.b - options.ineq.A*bnds(:,1);
end
% ===========================================

% ===========================================
% Context variables
% ===========================================
if isfield(options, 'cxtfun')
    h = options.cxtfun;
    h_len = length(h);
    if isa(h, 'function_handle')
        h = {h};
    end
elseif isfield(options, 'cxtlen')
    h_len = options.cxtlen;
else
    h_len = 0;
    h = {};
end

% My new option 'cxtmode' describes which approach SMGO uses to
% treat contextual optimization
% 1 - DIRECTLY using Mode 2, having the merit function decided by a mixture
%     of lower bounds and uncertainty
% 2 - Using Mode 1 and Mode 2, similar to classical SMGO
if isfield(options, 'cxtmode')
    cxt_mode = options.cxtmode;
else
    cxt_mode = CXTMODE_CLASSC;
end

% New option 'cxtrad' sets the ball radius around the current context, such that
% samples with contexts within this ball are "kind of in the same" context.
% For example, when deciding the best sample, all samples "in the same context"
% as the current one are collected, and the one with minimum objective and satisfying
% constraints is the "best sample". Furthermore, when the next context is within
% the same ball as the current one, it's treated as "kind of the same context",
% so we can do operations involving a trust region in the input variable space
if isfield(options, 'cxtrad')
    cxt_rad = options.cxtrad;
else
    cxt_rad = 0.1;
end
% ===========================================

% ===========================================
% Quasi-random points distribution parameters
% ===========================================
% TODO: sbl_seq should not just be the first sbl_size elements of sobolset
if isfield(options, 'sobolsize')
    sbl_size = options.sobolsize;
else
    sbl_size = 500;
end

if isfield(options, 'sobol')
    if ~options.sobol
        sbl_size = 0;
        sbl_seq  = [];
    end
end

if sbl_size
    try
        sbl_seq = sobolset(D-h_len);
        sbl_seq = sbl_seq(1:sbl_size, :)';
    catch
        sbl_size = 0;
        sbl_seq  = [];
    end
end

if sbl_size
    sbl_seq(:,2) = [];
    sbl_size     = sbl_size - 1;

    % Filtering Sobol-generated candidate points to eliminate 
    % those which are outside the linear inequalities
    if isfield(options,'ineq')
        sbl_rm = [];
        for sbl_i = 1:size(sbl_seq,2)
            if any(A_iq*sbl_seq(:,sbl_i)>b_iq,'all')
                sbl_rm = [sbl_rm sbl_i];
            end
        end
        sbl_seq(:,sbl_rm) = [];
        sbl_size          = sbl_size - length(sbl_rm);
    end
end

% ===========================================

% ===========================================
% Trust region parameters
% ===========================================
tr_size = 0.1;
if isfield(options, 'trustregion')
    if ~options.trustregion
        tr_size = 5.0;
    end
end
tr_coeff = 0.5;
% ===========================================

% ===========================================
% Noise bounds
% ===========================================
if isfield(options, 'noise')
    eps_col = options.noise;
else
    eps_col = NaN(1, 1+g_len);
end
% ===========================================

% ===========================================
% Initial samples
% ===========================================
if isfield(options, 'startX')
    X0      = options.startX;
    X0_rows = size(X0,1);
    if X0_rows ~= (D + 1 + g_len)
        fprintf("[ERROR] Incorrect data dimensions for options.startX. Please refer to documentation.\n");
        err_flag = true;
    end
    X0_real = real2normd(X0(1:D,:), bnds);
    if any(X0_real > 1.0, 'all') || any(X0_real < 0.0, 'all')
        fprintf("[WARNING] One or more points in options.startX are out of bounds. Will ignore such points.\n");        
        samp_rem_idx = (any(X0_real > 1.0, 1) | any(X0_real < 0.0, 1));
        X0(:, samp_rem_idx) = [];
    end
end
% ===========================================

% ===========================================
% Starting point
% ===========================================
if isfield(options, 'startpt')
    if strcmp(options.startpt, 'random')        
        x0 = rand(D, 1);
        if isfield(options,'ineq')
            rndgen_ctr = 0;
            while any(A_iq*x0 > b_iq,'all') && (rndgen_ctr < D*1000)
                x0 = rand(D, 1);
                rndgen_ctr = rndgen_ctr+1;
            end
            if A_iq*x0 <= b_iq
                disp('[INFO] Starting from a random start point.');
            else
                disp('[ERROR] Failed to generate random start point given the linear inequality. Are you sure that the resulting search space is not empty?');
                err_flag = true;
            end
        else
            disp('[INFO] Starting from a random start point.');
        end
    else
        x0 = real2normd(options.startpt, bnds);
        if isfield(options,'ineq')
            % TODO: if outside of the linear inequalities, do a projection
            if any(A_iq*x0 > b_iq,'all')
                disp('[ERROR] Specified starting point is outside the bounds! Have you checked if the start point is also within the linear inequalities (if present)?');
            end
            rndgen_ctr = 0;
            while any(A_iq*x0 > b_iq,'all') && (rndgen_ctr < D*1000)
                x0 = rand(D, 1);
                rndgen_ctr = rndgen_ctr+1;
            end
            if A_iq*x0 <= b_iq
                disp('[INFO] Starting instead from a randomly-generated start point.');
            else
                disp('[ERROR] Failed to generate random start point given the linear inequality. Are you sure that the resulting search space is not empty?');
                err_flag = true;
            end
        end
    end
else
    if exist('X0','var') && ~isempty(X0)
        x0 = [];
        disp('[INFO] No start point, but initial samples supplied. SMGO will immediately compute the next sampling point from initial samples.');
    else     
        x0 = rand(D, 1);
        if isfield(options,'ineq')
            rndgen_ctr = 0;
            while any(A_iq*x0 > b_iq,'all') && (rndgen_ctr < D*1000)
                x0 = rand(D, 1);
                rndgen_ctr = rndgen_ctr+1;
            end
            if A_iq*x0 <= b_iq
                disp('[INFO] Starting from a random start point.');
            else
                disp('[ERROR] Failed to generate random start point given the linear inequality. Are you sure that the resulting search space is not empty?');
                err_flag = true;
            end
        else
            disp('[INFO] Starting from a random start point.');
        end
    end
end
% If startpt is empty, then IF no initial samples, start from RANDOM
% If startpt is empty, and THERE ARE initial samples, start from those initial samples
% If startpt is NOT empty, then OK
% ===========================================

% ===========================================
% Feedback functions
% ===========================================
% Per-iteration feedback function
if isfield(options, 'iterfdbackfun')
    iterfdbackfun = options.iterfdbackfun;
else
    iterfdbackfun = {};
end

% Output feedback function
if isfield(options, 'outfdbackfun')
    outfdbackfun = options.outfdbackfun;
else
    outfdbackfun = {};
end

% ===========================================

% ===========================================
% SMGO-D advanced parameters
% ===========================================
if isfield(options, 'alpha')
    salpha = options.alpha;
else
    salpha = 0.005;
end

if isfield(options, 'beta')
    sbeta = options.beta;
else
    sbeta = 0.1;
end

if isfield(options, 'delta')
    sdelta = options.delta;
else
    sdelta = 0.2;
end

if isfield(options, 'phi')
    phi = options.phi;
else
    phi = 1e-6;
end

if isfield(options, 'B')
    B = options.B;
else
    B = 5;
end

% To enable filtering out of candidate points after iteration defined by 'filtercdpt_horiz'
% We do this to further speed up the algorithm, and to avoid too much memory usage
if isfield(options, 'filtercdpt')
    filtercdpt = options.filtercdpt;
else
    filtercdpt = true;
end

if isfield(options, 'filtercdpt_horiz')
    filtercdpt_horiz = options.filtercdpt_horiz;
else
    filtercdpt_horiz = 50*(D-h_len);
end

% ===========================================

% ============================================
% SMGO parameters for time-varying functions
% ============================================

% Set as 'true' if the hidden function is known as time-varying
if isfield(options, 'timevarying')
    timev = options.timevarying;
else
    timev = false;
end

% if isfield(options, 'initmodelfun')
%     timev_init = options.initmodelsamps;
% else
%     timev_init = 50*D;
% end
% 
% % The number of 'initial samples' that SMGO will get, for which the underlying function
% % is assumed to be in its 'initial state' (not varying with time yet).
% % 'timev_init' can be set to 0, and initial samples can be supplied through 'options.initX' instead
% if isfield(options, 'initmodelsamps')
%     timev_init = options.initmodelsamps;
% else
%     timev_init = 50*D;
% end

% The minimum 'age' of samples, after which they can be up for informativeness-based discarding
% NOTE: this is in terms of iteration number, not really the wall-clock time
if isfield(options, 'T_old')
    T_old = options.T_old;
else
    T_old = 100*D;
end

% The 'age' of samples, after which they are FORCE-discarded (because they are 'too old')
if isfield(options, 'T_oold')
    T_oold = options.T_oold;
else
    T_oold = 200*D;
end
% ============================================
