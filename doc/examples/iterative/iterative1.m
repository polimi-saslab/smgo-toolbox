%% =======================================
%  Iterative test #1: Styblinski-Tang function
%  [WARNING] Still a feature UNDER CONSTRUCTION
%  - we try to run an iterative version of SMGO
%    where we supply the initial samples, and
%    SMGO outputs the next point to sample as out.nxt_x.
%  - the out.nxt_x can be evaluated outside, and
%    iteratively added to the existing data set,
%    to be fed back to SMGO
%  =======================================

clear options; close all;

% Populating basic fields in options
options.objfun  = @(x)fn_styblinski(x);
options.bounds  = repmat([-5 5], 2, 1);
options.maxiter = 1; % For iterative implementation, we only need 1 SMGO iteration

% Generating 20 initial samples of the Styblinski-Tang function
% The array X0 should be of size (D + 1 + S) x N0 where
% - each column is a sample i
% - rows 1 to D are the location for sample i (D is the dimensionality)
% - row D + 1 is the objective sampled value
% - next S rows are the respective constraints sampled values
%   (S can be 0 for without black-box constraints, as in this example)
X0 = [];
for i = 1:50
    x0_n = normd2real(rand(2,1),options.bounds);
    X0 = [X0 [x0_n; fn_styblinski(x0_n)]];
end

% Putting the initial samples inside options.startX
options.startX = X0;

for cyc = 1:250
    % Run SMGO only for ONE iteration, just to get nxt_x
    out = smgo(options);

    % Now we have access to the next sampling point out.nxt_x, which we
    % will evaluate using an external function
    nxt_z = fn_styblinski(out.nxt_x);
    disp([out.nxt_x' nxt_z]);
    
    % Append our new sample to the previously-existing data set, and 
    % use it as our new initial samples options.startX
    options.startX = [options.startX [out.nxt_x; nxt_z]];
end

% Plot output
plot(options.startX(1,:), options.startX(2,:), '.'); grid on, title('Sampled points');