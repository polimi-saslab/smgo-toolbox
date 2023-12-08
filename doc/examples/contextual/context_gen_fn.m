%% =========================================
% Initializing generator of context sequence
% ==========================================

cxt_rand  = @(D)cxt_gen(D);
cxt_move  = @(D)cxt_gen_mv(D);
cxt_const = @(D)cxt_gen_const(D);

function cxt = cxt_gen(D)
    persistent rnd_gen 
    if isempty(rnd_gen)
        rnd_gen = rng(0, 'twister'); % Random generator 1
    end
    cxt = rand(D,1);
end

function cxt = cxt_gen_mv(D)
    persistent rnd_gen cxt_prev
    if isempty(rnd_gen)
        rnd_gen = rng(0, 'twister'); % Random generator 1
    end
    if isempty(cxt_prev)
        cxt_prev = rand(D,1);
    end
    cxt      = max(0.0, min(1.0, cxt_prev + 0.1*(rand(D,1)-0.5)));
    cxt_prev = cxt;
end

function cxt = cxt_gen_const(D)
    cxt = 0.5*ones(D,1);
end