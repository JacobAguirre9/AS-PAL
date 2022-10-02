%% SPDX-License-Identifier: MIT
% Copyright © 2021 Weiwei "William" Kong

% Solve a multivariate nonconvex quadratic programming problem constrained to the unit spectraplex intersected with an affine manifold.

% The function of interest is
%
%  f(x) :=  -xi / 2 * ||D * B * x|| ^ 2 + tau / 2 * ||A * x - b|| ^ 2
%
% with curvature pair (m, M).


% -------------------------------------------------------------------------
%% Global Variables
% -------------------------------------------------------------------------

% .........................................................................
% Create basic hparams.
base_hparam = struct();

ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;

qp_aipp_hparam = base_hparam;
qp_aipp_hparam.aipp_type = 'aipp';
qp_aipp_hparam.acg_steptype = 'variable';
qp_aipp_hparam.i_reset_prox_center = true;

rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';
rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;

iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;

% End basic hparams.
% .........................................................................

%%%%%%%%%%% Problem 1
Problem=1
% Create global hyperparams
N = 1000;
seed = 513;
dimM = 30;
dimN = 100;
M=10^1;
m=10^0;
density = 0.05;
global_tol = 1e-6;
time_limit = 14400;
% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_lin_cone_constr_02(N, M, m, seed, dimM, dimN, density);

% Create the Model object and specify the solver.
ncvx_lc_qp = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
ncvx_lc_qp.M = hparams.M;
ncvx_lc_qp.m = hparams.m;
ncvx_lc_qp.x0 = hparams.x0;
ncvx_lc_qp.K_constr = hparams.K_constr;

% Add linear constraints.
ncvx_lc_qp.constr_fn = @(x) hparams.constr_fn(x);
ncvx_lc_qp.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
ncvx_lc_qp.opt_type = 'relative';
ncvx_lc_qp.feas_type = 'relative';
ncvx_lc_qp.opt_tol = global_tol;
ncvx_lc_qp.feas_tol = global_tol;
ncvx_lc_qp.time_limit = time_limit;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {rqp_aipp_hparam, ialm_hparam};
name_arr = {'RQP', 'iALM'};
framework_arr = {@penalty3, @iALM3};
solver_arr = {@AIPP, @ECG};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);


RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);


ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);


filename='106SDPwk_June22.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');

%%%%%%%%%%% Problem 2
Problem=2
% Create global hyperparams
N = 1000;
seed = 164;
dimM = 30;
dimN = 100;
M=10^2;
m=10^0;
density = 0.05;
global_tol = 1e-6;
time_limit = 14400;
% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_lin_cone_constr_02(N, M, m, seed, dimM, dimN, density);

% Create the Model object and specify the solver.
ncvx_lc_qp = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
ncvx_lc_qp.M = hparams.M;
ncvx_lc_qp.m = hparams.m;
ncvx_lc_qp.x0 = hparams.x0;
ncvx_lc_qp.K_constr = hparams.K_constr;

% Add linear constraints.
ncvx_lc_qp.constr_fn = @(x) hparams.constr_fn(x);
ncvx_lc_qp.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
ncvx_lc_qp.opt_type = 'relative';
ncvx_lc_qp.feas_type = 'relative';
ncvx_lc_qp.opt_tol = global_tol;
ncvx_lc_qp.feas_tol = global_tol;
ncvx_lc_qp.time_limit = time_limit;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {rqp_aipp_hparam, ialm_hparam};
name_arr = {'RQP', 'iALM'};
framework_arr = {@penalty3, @iALM3};
solver_arr = {@AIPP, @ECG};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);


RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);


ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);


filename='106SDPwk_June22.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');

%%%%%%%%%%% Problem 3
Problem=3
% Create global hyperparams
N = 1000;
seed = 227;
dimM = 30;
dimN = 100;
M=10^3;
m=10^0;
density = 0.05;
global_tol = 1e-6;
time_limit = 14400;
% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_lin_cone_constr_02(N, M, m, seed, dimM, dimN, density);

% Create the Model object and specify the solver.
ncvx_lc_qp = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
ncvx_lc_qp.M = hparams.M;
ncvx_lc_qp.m = hparams.m;
ncvx_lc_qp.x0 = hparams.x0;
ncvx_lc_qp.K_constr = hparams.K_constr;

% Add linear constraints.
ncvx_lc_qp.constr_fn = @(x) hparams.constr_fn(x);
ncvx_lc_qp.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
ncvx_lc_qp.opt_type = 'relative';
ncvx_lc_qp.feas_type = 'relative';
ncvx_lc_qp.opt_tol = global_tol;
ncvx_lc_qp.feas_tol = global_tol;
ncvx_lc_qp.time_limit = time_limit;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {rqp_aipp_hparam, ialm_hparam};
name_arr = {'RQP', 'iALM'};
framework_arr = {@penalty3, @iALM3};
solver_arr = {@AIPP, @ECG};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);


RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);


ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);


filename='106SDPwk_June22.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');

%%%%%%%%%%% Problem 4
Problem=4
% Create global hyperparams
N = 1000;
seed = 622;
dimM = 30;
dimN = 100;
M=10^4;
m=10^0;
density = 0.05;
global_tol = 1e-6;
time_limit = 14400;
% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_lin_cone_constr_02(N, M, m, seed, dimM, dimN, density);

% Create the Model object and specify the solver.
ncvx_lc_qp = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
ncvx_lc_qp.M = hparams.M;
ncvx_lc_qp.m = hparams.m;
ncvx_lc_qp.x0 = hparams.x0;
ncvx_lc_qp.K_constr = hparams.K_constr;

% Add linear constraints.
ncvx_lc_qp.constr_fn = @(x) hparams.constr_fn(x);
ncvx_lc_qp.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
ncvx_lc_qp.opt_type = 'relative';
ncvx_lc_qp.feas_type = 'relative';
ncvx_lc_qp.opt_tol = global_tol;
ncvx_lc_qp.feas_tol = global_tol;
ncvx_lc_qp.time_limit = time_limit;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {rqp_aipp_hparam, ialm_hparam};
name_arr = {'RQP', 'iALM'};
framework_arr = {@penalty3, @iALM3};
solver_arr = {@AIPP, @ECG};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);


RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);


ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);


filename='106SDPwk_June22.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');


%%%%%%%%%%% Problem 5
Problem=5
% Create global hyperparams
N = 1000;
seed = 326;
dimM = 30;
dimN = 100;
M=10^5;
m=10^0;
density = 0.05;
global_tol = 1e-6;
time_limit = 14400;
% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_lin_cone_constr_02(N, M, m, seed, dimM, dimN, density);

% Create the Model object and specify the solver.
ncvx_lc_qp = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
ncvx_lc_qp.M = hparams.M;
ncvx_lc_qp.m = hparams.m;
ncvx_lc_qp.x0 = hparams.x0;
ncvx_lc_qp.K_constr = hparams.K_constr;

% Add linear constraints.
ncvx_lc_qp.constr_fn = @(x) hparams.constr_fn(x);
ncvx_lc_qp.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
ncvx_lc_qp.opt_type = 'relative';
ncvx_lc_qp.feas_type = 'relative';
ncvx_lc_qp.opt_tol = global_tol;
ncvx_lc_qp.feas_tol = global_tol;
ncvx_lc_qp.time_limit = time_limit;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {rqp_aipp_hparam, ialm_hparam};
name_arr = {'RQP', 'iALM'};
framework_arr = {@penalty3, @iALM3};
solver_arr = {@AIPP, @ECG};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);


RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);


ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);


filename='106SDPwk_June22.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');


%%%%%%%%%%% Problem 6
Problem=6
% Create global hyperparams
N = 1000;
seed = 648;
dimM = 30;
dimN = 100;
M=10^2;
m=10^1;
density = 0.05;
global_tol = 1e-6;
time_limit = 14400;
% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_lin_cone_constr_02(N, M, m, seed, dimM, dimN, density);

% Create the Model object and specify the solver.
ncvx_lc_qp = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
ncvx_lc_qp.M = hparams.M;
ncvx_lc_qp.m = hparams.m;
ncvx_lc_qp.x0 = hparams.x0;
ncvx_lc_qp.K_constr = hparams.K_constr;

% Add linear constraints.
ncvx_lc_qp.constr_fn = @(x) hparams.constr_fn(x);
ncvx_lc_qp.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
ncvx_lc_qp.opt_type = 'relative';
ncvx_lc_qp.feas_type = 'relative';
ncvx_lc_qp.opt_tol = global_tol;
ncvx_lc_qp.feas_tol = global_tol;
ncvx_lc_qp.time_limit = time_limit;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {rqp_aipp_hparam, ialm_hparam};
name_arr = {'RQP', 'iALM'};
framework_arr = {@penalty3, @iALM3};
solver_arr = {@AIPP, @ECG};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);


RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);

filename='106SDPwk_June22.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');

%%%%%%%%%%% Problem 7
Problem=7
% Create global hyperparams
N = 1000;
seed = 734;
dimM = 30;
dimN = 100;
M=10^3;
m=10^1;
density = 0.05;
global_tol = 1e-6;
time_limit = 14400;
% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_lin_cone_constr_02(N, M, m, seed, dimM, dimN, density);

% Create the Model object and specify the solver.
ncvx_lc_qp = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
ncvx_lc_qp.M = hparams.M;
ncvx_lc_qp.m = hparams.m;
ncvx_lc_qp.x0 = hparams.x0;
ncvx_lc_qp.K_constr = hparams.K_constr;

% Add linear constraints.
ncvx_lc_qp.constr_fn = @(x) hparams.constr_fn(x);
ncvx_lc_qp.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
ncvx_lc_qp.opt_type = 'relative';
ncvx_lc_qp.feas_type = 'relative';
ncvx_lc_qp.opt_tol = global_tol;
ncvx_lc_qp.feas_tol = global_tol;
ncvx_lc_qp.time_limit = time_limit;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {rqp_aipp_hparam, ialm_hparam};
name_arr = {'RQP', 'iALM'};
framework_arr = {@penalty3, @iALM3};
solver_arr = {@AIPP, @ECG};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);


RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);


ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);


filename='106SDPwk_June22.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');

%%%%%%%%%%% Problem 8
Problem=8
% Create global hyperparams
N = 1000;
seed = 303;
dimM = 30;
dimN = 100;
M=10^4;
m=10^1;
density = 0.05;
global_tol = 1e-6;
time_limit = 14400;
% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_lin_cone_constr_02(N, M, m, seed, dimM, dimN, density);

% Create the Model object and specify the solver.
ncvx_lc_qp = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
ncvx_lc_qp.M = hparams.M;
ncvx_lc_qp.m = hparams.m;
ncvx_lc_qp.x0 = hparams.x0;
ncvx_lc_qp.K_constr = hparams.K_constr;

% Add linear constraints.
ncvx_lc_qp.constr_fn = @(x) hparams.constr_fn(x);
ncvx_lc_qp.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
ncvx_lc_qp.opt_type = 'relative';
ncvx_lc_qp.feas_type = 'relative';
ncvx_lc_qp.opt_tol = global_tol;
ncvx_lc_qp.feas_tol = global_tol;
ncvx_lc_qp.time_limit = time_limit;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {rqp_aipp_hparam, ialm_hparam};
name_arr = {'RQP', 'iALM'};
framework_arr = {@penalty3, @iALM3};
solver_arr = {@AIPP, @ECG};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);


RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);


ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);


filename='106SDPwk_June22.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');


%%%%%%%%%%% Problem 9
Problem=9
% Create global hyperparams
N = 1000;
seed = 871;
dimM = 30;
dimN = 100;
M=10^5;
m=10^1;
density = 0.05;
global_tol = 1e-6;
time_limit = 14400;
% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_lin_cone_constr_02(N, M, m, seed, dimM, dimN, density);

% Create the Model object and specify the solver.
ncvx_lc_qp = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
ncvx_lc_qp.M = hparams.M;
ncvx_lc_qp.m = hparams.m;
ncvx_lc_qp.x0 = hparams.x0;
ncvx_lc_qp.K_constr = hparams.K_constr;

% Add linear constraints.
ncvx_lc_qp.constr_fn = @(x) hparams.constr_fn(x);
ncvx_lc_qp.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
ncvx_lc_qp.opt_type = 'relative';
ncvx_lc_qp.feas_type = 'relative';
ncvx_lc_qp.opt_tol = global_tol;
ncvx_lc_qp.feas_tol = global_tol;
ncvx_lc_qp.time_limit = time_limit;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {rqp_aipp_hparam, ialm_hparam};
name_arr = {'RQP', 'iALM'};
framework_arr = {@penalty3, @iALM3};
solver_arr = {@AIPP, @ECG};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);


RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);


ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);


filename='106SDPwk_June22.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');

%%%%%%%%%%% Problem 10
Problem=10
% Create global hyperparams
N = 1000;
seed = 535;
dimM = 30;
dimN = 100;
M=10^6;
m=10^1;
density = 0.05;
global_tol = 1e-6;
time_limit = 14400;
% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_lin_cone_constr_02(N, M, m, seed, dimM, dimN, density);

% Create the Model object and specify the solver.
ncvx_lc_qp = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
ncvx_lc_qp.M = hparams.M;
ncvx_lc_qp.m = hparams.m;
ncvx_lc_qp.x0 = hparams.x0;
ncvx_lc_qp.K_constr = hparams.K_constr;

% Add linear constraints.
ncvx_lc_qp.constr_fn = @(x) hparams.constr_fn(x);
ncvx_lc_qp.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
ncvx_lc_qp.opt_type = 'relative';
ncvx_lc_qp.feas_type = 'relative';
ncvx_lc_qp.opt_tol = global_tol;
ncvx_lc_qp.feas_tol = global_tol;
ncvx_lc_qp.time_limit = time_limit;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {rqp_aipp_hparam, ialm_hparam};
name_arr = {'RQP', 'iALM'};
framework_arr = {@penalty3, @iALM3};
solver_arr = {@AIPP, @ECG};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);


RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);


ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);


filename='106SDPwk_June22.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');


%%%%%%%%%%% Problem 11
Problem=11
% Create global hyperparams
N = 1000;
seed = 421;
dimM = 30;
dimN = 100;
M=10^4;
m=10^2;
density = 0.05;
global_tol = 1e-6;
time_limit = 14400;
% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_lin_cone_constr_02(N, M, m, seed, dimM, dimN, density);

% Create the Model object and specify the solver.
ncvx_lc_qp = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
ncvx_lc_qp.M = hparams.M;
ncvx_lc_qp.m = hparams.m;
ncvx_lc_qp.x0 = hparams.x0;
ncvx_lc_qp.K_constr = hparams.K_constr;

% Add linear constraints.
ncvx_lc_qp.constr_fn = @(x) hparams.constr_fn(x);
ncvx_lc_qp.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
ncvx_lc_qp.opt_type = 'relative';
ncvx_lc_qp.feas_type = 'relative';
ncvx_lc_qp.opt_tol = global_tol;
ncvx_lc_qp.feas_tol = global_tol;
ncvx_lc_qp.time_limit = time_limit;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {rqp_aipp_hparam, ialm_hparam};
name_arr = {'RQP', 'iALM'};
framework_arr = {@penalty3, @iALM3};
solver_arr = {@AIPP, @ECG};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);


RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);


ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);


filename='106SDPwk_June22.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');

%%%%%%%%%%% Problem 12
Problem=12
% Create global hyperparams
N = 1000;
seed = 310;
dimM = 30;
dimN = 100;
M=10^6;
m=10^2;
density = 0.05;
global_tol = 1e-6;
time_limit = 14400;
% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_lin_cone_constr_02(N, M, m, seed, dimM, dimN, density);

% Create the Model object and specify the solver.
ncvx_lc_qp = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
ncvx_lc_qp.M = hparams.M;
ncvx_lc_qp.m = hparams.m;
ncvx_lc_qp.x0 = hparams.x0;
ncvx_lc_qp.K_constr = hparams.K_constr;

% Add linear constraints.
ncvx_lc_qp.constr_fn = @(x) hparams.constr_fn(x);
ncvx_lc_qp.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
ncvx_lc_qp.opt_type = 'relative';
ncvx_lc_qp.feas_type = 'relative';
ncvx_lc_qp.opt_tol = global_tol;
ncvx_lc_qp.feas_tol = global_tol;
ncvx_lc_qp.time_limit = time_limit;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {rqp_aipp_hparam, ialm_hparam};
name_arr = {'RQP', 'iALM'};
framework_arr = {@penalty3, @iALM3};
solver_arr = {@AIPP, @ECG};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);


RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);


ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);


filename='106SDPwk_June22.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');

%%%%%%%%%%% Problem 13
Problem=13
% Create global hyperparams
N = 1000;
seed = 779;
dimM = 30;
dimN = 100;
M=10^7;
m=10^2;
density = 0.05;
global_tol = 1e-6;
time_limit = 14400;
% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_lin_cone_constr_02(N, M, m, seed, dimM, dimN, density);

% Create the Model object and specify the solver.
ncvx_lc_qp = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
ncvx_lc_qp.M = hparams.M;
ncvx_lc_qp.m = hparams.m;
ncvx_lc_qp.x0 = hparams.x0;
ncvx_lc_qp.K_constr = hparams.K_constr;

% Add linear constraints.
ncvx_lc_qp.constr_fn = @(x) hparams.constr_fn(x);
ncvx_lc_qp.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
ncvx_lc_qp.opt_type = 'relative';
ncvx_lc_qp.feas_type = 'relative';
ncvx_lc_qp.opt_tol = global_tol;
ncvx_lc_qp.feas_tol = global_tol;
ncvx_lc_qp.time_limit = time_limit;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {rqp_aipp_hparam, ialm_hparam};
name_arr = {'RQP', 'iALM'};
framework_arr = {@penalty3, @iALM3};
solver_arr = {@AIPP, @ECG};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);


RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);


ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);


filename='106SDPwk_June22.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');


%%%%%%%%%%% Problem 14
Problem=14
% Create global hyperparams
N = 1000;
seed = 813;
dimM = 30;
dimN = 100;
M=10^4;
m=10^3;
density = 0.05;
global_tol = 1e-6;
time_limit = 14400;
% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_lin_cone_constr_02(N, M, m, seed, dimM, dimN, density);

% Create the Model object and specify the solver.
ncvx_lc_qp = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
ncvx_lc_qp.M = hparams.M;
ncvx_lc_qp.m = hparams.m;
ncvx_lc_qp.x0 = hparams.x0;
ncvx_lc_qp.K_constr = hparams.K_constr;

% Add linear constraints.
ncvx_lc_qp.constr_fn = @(x) hparams.constr_fn(x);
ncvx_lc_qp.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
ncvx_lc_qp.opt_type = 'relative';
ncvx_lc_qp.feas_type = 'relative';
ncvx_lc_qp.opt_tol = global_tol;
ncvx_lc_qp.feas_tol = global_tol;
ncvx_lc_qp.time_limit = time_limit;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {rqp_aipp_hparam, ialm_hparam};
name_arr = {'RQP', 'iALM'};
framework_arr = {@penalty3, @iALM3};
solver_arr = {@AIPP, @ECG};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);


RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);


ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);


filename='106SDPwk_June22.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');


%%%%%%%%%%% Problem 15
Problem=15
% Create global hyperparams
N = 1000;
seed = 544;
dimM = 30;
dimN = 100;
M=10^5;
m=10^3;
density = 0.05;
global_tol = 1e-6;
time_limit = 14400;
% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_lin_cone_constr_02(N, M, m, seed, dimM, dimN, density);

% Create the Model object and specify the solver.
ncvx_lc_qp = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
ncvx_lc_qp.M = hparams.M;
ncvx_lc_qp.m = hparams.m;
ncvx_lc_qp.x0 = hparams.x0;
ncvx_lc_qp.K_constr = hparams.K_constr;

% Add linear constraints.
ncvx_lc_qp.constr_fn = @(x) hparams.constr_fn(x);
ncvx_lc_qp.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
ncvx_lc_qp.opt_type = 'relative';
ncvx_lc_qp.feas_type = 'relative';
ncvx_lc_qp.opt_tol = global_tol;
ncvx_lc_qp.feas_tol = global_tol;
ncvx_lc_qp.time_limit = time_limit;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {rqp_aipp_hparam, ialm_hparam};
name_arr = {'RQP', 'iALM'};
framework_arr = {@penalty3, @iALM3};
solver_arr = {@AIPP, @ECG};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);


RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);


ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);


filename='106SDPwk_June22.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');





%%%%%%%%%%% Problem 16
Problem=16
% Create global hyperparams
N = 1000;
seed = 881;
dimM = 30;
dimN = 100;
M=10^7;
m=10^3;
density = 0.05;
global_tol = 1e-6;
time_limit = 14400;
% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_lin_cone_constr_02(N, M, m, seed, dimM, dimN, density);

% Create the Model object and specify the solver.
ncvx_lc_qp = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
ncvx_lc_qp.M = hparams.M;
ncvx_lc_qp.m = hparams.m;
ncvx_lc_qp.x0 = hparams.x0;
ncvx_lc_qp.K_constr = hparams.K_constr;

% Add linear constraints.
ncvx_lc_qp.constr_fn = @(x) hparams.constr_fn(x);
ncvx_lc_qp.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
ncvx_lc_qp.opt_type = 'relative';
ncvx_lc_qp.feas_type = 'relative';
ncvx_lc_qp.opt_tol = global_tol;
ncvx_lc_qp.feas_tol = global_tol;
ncvx_lc_qp.time_limit = time_limit;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {rqp_aipp_hparam, ialm_hparam};
name_arr = {'RQP', 'iALM'};
framework_arr = {@penalty3, @iALM3};
solver_arr = {@AIPP, @ECG};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);


RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);


ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);


filename='106SDPwk_June22.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');

%%%%%%%%%%% Problem 17
Problem=17
% Create global hyperparams
N = 1000;
seed = 121;
dimM = 30;
dimN = 100;
M=10^8;
m=10^3;
density = 0.05;
global_tol = 1e-6;
time_limit = 14400;
% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_lin_cone_constr_02(N, M, m, seed, dimM, dimN, density);

% Create the Model object and specify the solver.
ncvx_lc_qp = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
ncvx_lc_qp.M = hparams.M;
ncvx_lc_qp.m = hparams.m;
ncvx_lc_qp.x0 = hparams.x0;
ncvx_lc_qp.K_constr = hparams.K_constr;

% Add linear constraints.
ncvx_lc_qp.constr_fn = @(x) hparams.constr_fn(x);
ncvx_lc_qp.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
ncvx_lc_qp.opt_type = 'relative';
ncvx_lc_qp.feas_type = 'relative';
ncvx_lc_qp.opt_tol = global_tol;
ncvx_lc_qp.feas_tol = global_tol;
ncvx_lc_qp.time_limit = time_limit;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {rqp_aipp_hparam, ialm_hparam};
name_arr = {'RQP', 'iALM'};
framework_arr = {@penalty3, @iALM3};
solver_arr = {@AIPP, @ECG};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);


RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);


ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);


filename='106SDPwk_June22.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');


%%%%%%%%%%% Problem 18
Problem=18
% Create global hyperparams
N = 1000;
seed = 217;
dimM = 30;
dimN = 100;
M=10^5;
m=10^4;
density = 0.05;
global_tol = 1e-6;
time_limit = 14400;
% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_lin_cone_constr_02(N, M, m, seed, dimM, dimN, density);

% Create the Model object and specify the solver.
ncvx_lc_qp = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
ncvx_lc_qp.M = hparams.M;
ncvx_lc_qp.m = hparams.m;
ncvx_lc_qp.x0 = hparams.x0;
ncvx_lc_qp.K_constr = hparams.K_constr;

% Add linear constraints.
ncvx_lc_qp.constr_fn = @(x) hparams.constr_fn(x);
ncvx_lc_qp.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
ncvx_lc_qp.opt_type = 'relative';
ncvx_lc_qp.feas_type = 'relative';
ncvx_lc_qp.opt_tol = global_tol;
ncvx_lc_qp.feas_tol = global_tol;
ncvx_lc_qp.time_limit = time_limit;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {rqp_aipp_hparam, ialm_hparam};
name_arr = {'RQP', 'iALM'};
framework_arr = {@penalty3, @iALM3};
solver_arr = {@AIPP, @ECG};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);


RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);


ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);


filename='106SDPwk_June22.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');

%%%%%%%%%%% Problem 19
Problem=19
% Create global hyperparams
N = 1000;
seed = 218;
dimM = 30;
dimN = 100;
M=10^6;
m=10^4;
density = 0.05;
global_tol = 1e-6;
time_limit = 14400;
% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_lin_cone_constr_02(N, M, m, seed, dimM, dimN, density);

% Create the Model object and specify the solver.
ncvx_lc_qp = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
ncvx_lc_qp.M = hparams.M;
ncvx_lc_qp.m = hparams.m;
ncvx_lc_qp.x0 = hparams.x0;
ncvx_lc_qp.K_constr = hparams.K_constr;

% Add linear constraints.
ncvx_lc_qp.constr_fn = @(x) hparams.constr_fn(x);
ncvx_lc_qp.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
ncvx_lc_qp.opt_type = 'relative';
ncvx_lc_qp.feas_type = 'relative';
ncvx_lc_qp.opt_tol = global_tol;
ncvx_lc_qp.feas_tol = global_tol;
ncvx_lc_qp.time_limit = time_limit;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {rqp_aipp_hparam, ialm_hparam};
name_arr = {'RQP', 'iALM'};
framework_arr = {@penalty3, @iALM3};
solver_arr = {@AIPP, @ECG};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);


RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);


ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);


filename='106SDPwk_June22.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');

%%%%%%%%%%% Problem 20
Problem=20
% Create global hyperparams
N = 1000;
seed = 357;
dimM = 30;
dimN = 100;
M=10^8;
m=10^4;
density = 0.05;
global_tol = 1e-6;
time_limit = 14400;
% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_lin_cone_constr_02(N, M, m, seed, dimM, dimN, density);

% Create the Model object and specify the solver.
ncvx_lc_qp = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
ncvx_lc_qp.M = hparams.M;
ncvx_lc_qp.m = hparams.m;
ncvx_lc_qp.x0 = hparams.x0;
ncvx_lc_qp.K_constr = hparams.K_constr;

% Add linear constraints.
ncvx_lc_qp.constr_fn = @(x) hparams.constr_fn(x);
ncvx_lc_qp.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
ncvx_lc_qp.opt_type = 'relative';
ncvx_lc_qp.feas_type = 'relative';
ncvx_lc_qp.opt_tol = global_tol;
ncvx_lc_qp.feas_tol = global_tol;
ncvx_lc_qp.time_limit = time_limit;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {rqp_aipp_hparam, ialm_hparam};
name_arr = {'RQP', 'iALM'};
framework_arr = {@penalty3, @iALM3};
solver_arr = {@AIPP, @ECG};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);


RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);


ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);


filename='106SDPwk_June22.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');
























































































































