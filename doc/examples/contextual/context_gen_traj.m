% D_c        = 2;   % Dimensionality of the context
max_iter   = 1000; % Number of iteration, for which we generate the context
max_trials = 25;  % Number of optimization runs ("trials")
p_jump     = 0.1; % Probability at every iteration of generating a random ("jumping") context
                  % If not triggered, the context will just move at a random direction

for D_c = 1:10
    cxt_traj_col = [];
    for j = 1:max_trials
        cxt_traj = zeros(D_c, max_iter);
        cxt_traj(:,1) = rand(D_c,1);
        
        for i = 2:max_iter
            if rand(1) <= p_jump
                cxt_traj(:,i) = rand(D_c,1);
            else
                cxt_traj(:,i) = min(1, max(0, cxt_traj(:,i-1) + 0.05*(rand(D_c,1)-0.5)));
            end
        end
    
        cxt_traj_col = [cxt_traj_col; cxt_traj];
    end
    
    writematrix(cxt_traj_col, sprintf("cxt_%dD.csv",D_c));
end
