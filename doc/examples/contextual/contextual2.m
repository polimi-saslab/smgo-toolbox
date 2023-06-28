%% ============================================
% Contextual test #2: Styblinski-Tang function
% =============================================
% - 1st input is the (controllable) input variable
% - 2nd input is the (uncontrollable but observable) context variable
% - we will have 50 initial samples
%  ============================================

clear options; close all;

init_context_gen;

% Populating basic fields in options
options.objfun  = @(x)fn_styblinski(x);
options.bounds  = repmat([-5 5], 2, 1);
options.cxtfun  = @(x)normd2real(cxt_rand(1),[-5 5]);
options.maxiter = 500;
options.beta    = 0.1;
options.alpha   = 0.0;
options.filtercdpt = true;

X0 = [];
for i = 1:50
    x0_n = normd2real(rand(2,1),options.bounds);
    X0 = [X0 [x0_n; fn_styblinski(x0_n)]];
end

% Putting the initial samples inside options.startX
options.startX = X0;

% Run SMGO
out             = smgo(options);

% Output animation
animate_result;
