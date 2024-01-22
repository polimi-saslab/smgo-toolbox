%% ================================
%  Calculating noise
%  ================================

% Collecting which points are within radius 'rad' from each other
% Increment 'rad' when there are no points inside the ball of any other
for rad = 0.1:0.1:1.0
    ball_rad = rad * diam/2;
    ball_idx = zeros(X_n_len);
    for k = 1:X_n_len
        ball_idx(k, :) = vecnorm(X_n(1:D, :) - repmat(X_n(1:D, k), 1, X_n_len)) < ball_rad;
        ball_idx(k, k) = 0;
    end
    if sum(ball_idx,'all') > 0
        break;
    end
end
ball_idx = logical(ball_idx);
if iter == 3
    iter;
end
% Collecting the noise values
if isnan(eps_col(1))
    % If first entry of eps_col (supplied as options.noise) is NaN, let SMGO estimate the noise bound
    ns_lkup = NaN(1, X_n_len);
    for j = 1:X_n_len
        if sum(ball_idx(j, :)) > 0
            ns_lkup(j) = max(abs(X_n(XN_FVAL, ball_idx(j, :)) - X_n(XN_FVAL, j)));
        end
    end
    feps = sum(ns_lkup(ns_lkup > 0)) / (2 * sum(ball_idx,'all'));
else
    feps = eps_col(1);
end
for g_i = 1:g_len
    if isnan(eps_col(g_i + 1))
        % If entry of eps_col (supplied as options.noise) is NaN, let SMGO estimate the noise bound
        XN_GVAL_I = XN_GVAL + g_i - 1;            
        ns_lkup = NaN(1, X_n_len);
        for j = 1:X_n_len
            if sum(ball_idx(j, :)) > 0
                ns_lkup(j) = max(abs(X_n(XN_GVAL_I, ball_idx(j, :)) - X_n(XN_GVAL_I, j)));
            end            
        end
        geps(g_i) = sum(ns_lkup(ns_lkup > 0)) / (2 * sum(ball_idx,'all'));
    else
        geps(g_i) = eps_col(g_i + 1);
    end
end 