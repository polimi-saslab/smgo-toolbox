%% ============================================
% Contextual test #1: Styblinski-Tang function
% =============================================
% - 1st input is the (controllable) input variable
% - 2nd input is the (uncontrollable but observable) context variable
%  ============================================

clear options; close all;

init_context_gen;

% Populating basic fields in options
options.objfun  = @(x)fn_styblinski(x);
options.bounds  = repmat([-5 5], 2, 1);
options.cxtfun  = @(x)normd2real(cxt_rand(1),[-5 5]);
options.maxiter = 500;
options.beta    = 0.1;
options.alpha   = 0.005;
options.filtercdpt = true;
options.cxtmode = 2;

% Run SMGO
out             = smgo(options);

% Output animation
animate_result;