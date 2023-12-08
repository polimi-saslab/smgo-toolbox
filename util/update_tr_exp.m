if ~h_len && iter > 1
    if (mode_prev == MODE_EXPLOIT) && opt_z_upd
        if z_n < exploit_thr
            % If the actual sample is actually less than expected improvement, 
            % then enlarge the trust region hyperbox (decrease TR exponent)
            tr_exp = max(0, tr_exp - 1);
        end
    else
        % Shrink the trust region hyperbox
        if tr_exp < 10
            tr_exp = tr_exp + 1;
        elseif mode_prev == MODE_EXPLORE
            tr_exp = tr_exp_0;
        end
    end    
end
tr_bnds = ones(D, 1) * [-0.5 0.5] * tr_size * (tr_coeff ^ tr_exp) + opt_x * ones(1, 2);
if h_len
    tr_bnds(1:h_len,:) = repmat(h_n,1,2);
end
tr_bnds = max(min(tr_bnds, 1.0), 0.0);