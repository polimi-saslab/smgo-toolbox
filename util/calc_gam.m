%% ========================================
%  Generate/update Lipschitz constants
%  ========================================

fgam_prev = fgam; 
fgam_db = zeros(X_n_len, 1);
for i = 1:X_n_len
    elg_idx = abs(X_n(XN_FVAL, :) - X_n(XN_FVAL, i)) > 2*feps;            
    if sum(elg_idx) > 0
        fgam_db(i) = max((abs(X_n(XN_FVAL, elg_idx) - X_n(XN_FVAL, i)) - 2*feps) ./ ...
                          vecnorm(X_n(1:D, elg_idx) - repmat(X_n(1:D, i), 1, sum(elg_idx))));
    else 
        fgam_db(i) = 1e-6;
    end
end
fgam = max(fgam_db);

ggam_prev = ggam;
for g_i = 1:g_len
    ggam_db = zeros(X_n_len, 1);
    XN_GVAL_I = XN_GVAL + g_i - 1;
    for i = 1:X_n_len
        elg_idx = abs(X_n(XN_GVAL_I, :) - X_n(XN_GVAL_I, i)) > 2*geps(g_i);
        if sum(elg_idx) > 0
            ggam_db(i) = max((abs(X_n(XN_GVAL_I, elg_idx) - X_n(XN_GVAL_I, i)) - 2*geps(g_i)) ./ ...
                              vecnorm(X_n(1:D, elg_idx) - repmat(X_n(1:D, i), 1, sum(elg_idx))));
        else 
            ggam_db(i) = 1e-6;
        end
    end
    ggam(g_i) = max(ggam_db);
end