%% ====================================
% Row definitions for samples set (X_n)
% =====================================

XN_FVAL  = D + 1; % objective value
XN_GVAL  = D + 2; % constraints values (if present)
XN_ID    = D + 2 + g_len;
XN_TSTMP = D + 3 + g_len;
xn_id        = 1;

%% ==================================================
% Sobol sequence-based exploitation array declaration
% ===================================================
SBL_FUB    = D + 1;
SBL_FLB    = D + 2;
SBL_G_INFO = D + 3;
SBL_GUB    = 0;
SBL_GLB    = 1;

%% =====================================================================
% Candidate points database definitions (for iteratively-updated bounds)
% ======================================================================

if ~h_len
    % ==================================================
    % Candidate points database db_cdpt: array structure
    % ==================================================
    % - rows 1 to D: the cdpt location
    % - row D+1: the value of upper bound cone vertex at that cdpt
    % - row D+2: the value of upper bound cone height at that cdpt
    % - row D+3: the value of lower bound cone vertex at that cdpt
    % - row D+4: the value of lower bound cone height at that cdpt
    % - row D+5: the value of uncertainty (lambda) at that cdpt
    CDPT_F_UB_VTX = D + 1;
    CDPT_F_UB_HGT = D + 2;
    CDPT_F_LB_VTX = D + 3;
    CDPT_F_LB_HGT = D + 4;
    CDPT_F_LAMBDA = D + 5;
    CDPT_F_ROWS   = 5;
    
    % the succeeding rows of db_cdpt describe similar cone-related information, but for each
    %   constraint. As all constraints are black-box too, we estimate the bounds
    %   and the uncertainty (Pi_s)
    CDPT_G_INFO   = D + 6;
    CDPT_G_UB_VTX = 0;
    CDPT_G_UB_HGT = 1;
    CDPT_G_LB_VTX = 2;
    CDPT_G_LB_HGT = 3;
    CDPT_G_LAMBDA = 4;
    CDPT_G_EST    = 5;    
    CDPT_G_ROWS   = 6;
    
    % The last two rows describe the remoteness measure and the candidate point age
    % end-1: Remoteness measure of each candidate point
    % end:   Age of candidate point (measured as the difference between iteration 
    %        of its generation to the current iteration)
    % ===============================================
    
    % ===================================================================
    % Row definitions for cone-generator arrays db_cdpt_cn and db_samp_cn
    % ===================================================================
    CDPTCN_UB = 1;
    CDPTCN_LB = 2;
    
    % ==================================================
    % Sampled points database db_samp: array structure
    % ==================================================
    % - rows 1 to D: the cdpoint coordinate
    % - row D+1: the value of upper bound cone vertex at that cdpt
    % - row D+2: the value of upper bound cone height at that cdpt
    % - row D+3: the value of lower bound cone vertex at that cdpt
    % - row D+4: the value of lower bound cone height at that cdpt
    % - row D+5: the value of uncertainty (lambda) at that cdpt
    SAMP_F_UB_VTX = 1;
    SAMP_F_UB_HGT = 2;
    SAMP_F_LB_VTX = 3;
    SAMP_F_LB_HGT = 4;
    
    % The succeeding rows of Mode psi describe similar information, but for each
    %   constraint. As all constraints are black-box too, we estimate the bounds
    %   and the uncertainty (Pi_s)
    SAMP_G_INFO   = 5;
    SAMP_G_UB_VTX = 0;
    SAMP_G_UB_HGT = 1;
    SAMP_G_LB_VTX = 2;
    SAMP_G_LB_HGT = 3;
    % =================================================
end

%% =====================================================================
% Candidate points database definitions (for exactly-calculated bounds)
% ======================================================================
if h_len
    % ==================================================
    % Candidate points database db_cdpt: array structure
    % ==================================================
    CDPT_F_UB     = D + 1;
    CDPT_F_LB     = D + 2;
    CDPT_F_LAMBDA = D + 3;
    CDPT_F_ROWS   = 3;

    CDPT_G_INFO   = D + 3;
    CDPT_G_UB     = 0;
    CDPT_G_LB     = 1;
    CDPT_G_LAMBDA = 2;
    CDPT_G_EST    = 3;
    CDPT_G_ROWS   = 4;
    
    % The last two rows describe the remoteness measure and the candidate point age
    % end-1: Remoteness measure of each candidate point
    % end:   Age of candidate point (measured as the difference between iteration 
    %        of its generation to the current iteration)
    % ===============================================
end
