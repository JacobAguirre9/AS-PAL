%% SPDX-License-Identifier: MIT
% Copyright © 2021 Weiwei "William" Kong

% Solve a multivariate nonconvex quadratic programming problem constrained to the unit simplex intersected with an affine manifold.

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
format long
base_hparam = struct();

ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;


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

%%%%%% Problem 1
% Create global hyperparams
N = 1000;
seed = 118;
dimM = 20;
dimN = 1000;
M=10^1;
m=10^0;
global_tol = 1e-4;
time_limit = 10800;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

[oracle, hparams] = test_fn_lin_cone_constr_01(N, M, m, seed, dimM, dimN);

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
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL','RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);


FVAL(1,1)=summary_tables.fval(1,1);
FVAL(1,2)=summary_tables.fval(1,2);
FVAL(1,3)=summary_tables.fval(1,3);

RUN(1,1)=summary_tables.runtime(1,1);
RUN(1,2)=summary_tables.runtime(1,2);
RUN(1,3)=summary_tables.runtime(1,3);

ACGiter(1,1)=summary_tables.iter(1,1);
ACGiter(1,2)=summary_tables.iter(1,2);
ACGiter(1,3)=summary_tables.iter(1,3);


filename='104Vector14instances.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');








%%%%%% Problem 2
% Create global hyperparams
N = 1000;
seed = 661;
dimM = 20;
dimN = 1000;
M=10^2;
m=10^0;
global_tol = 1e-4;
time_limit = 10800;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

[oracle, hparams] = test_fn_lin_cone_constr_01(N, M, m, seed, dimM, dimN);

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
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL','RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);


FVAL(2,1)=summary_tables.fval(1,1);
FVAL(2,2)=summary_tables.fval(1,2);
FVAL(2,3)=summary_tables.fval(1,3);

RUN(2,1)=summary_tables.runtime(1,1);
RUN(2,2)=summary_tables.runtime(1,2);
RUN(2,3)=summary_tables.runtime(1,3);

ACGiter(2,1)=summary_tables.iter(1,1);
ACGiter(2,2)=summary_tables.iter(1,2);
ACGiter(2,3)=summary_tables.iter(1,3);

filename='104Vector14instances.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');


%%%%%%%%% Problem 3
% Create global hyperparams
N = 1000;
seed = 815;
dimM = 20;
dimN = 1000;
M=10^3;
m=10^0;
global_tol = 1e-4;
time_limit = 10800;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

[oracle, hparams] = test_fn_lin_cone_constr_01(N, M, m, seed, dimM, dimN);

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
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL','RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);


FVAL(3,1)=summary_tables.fval(1,1);
FVAL(3,2)=summary_tables.fval(1,2);
FVAL(3,3)=summary_tables.fval(1,3);

RUN(3,1)=summary_tables.runtime(1,1);
RUN(3,2)=summary_tables.runtime(1,2);
RUN(3,3)=summary_tables.runtime(1,3);

ACGiter(3,1)=summary_tables.iter(1,1);
ACGiter(3,2)=summary_tables.iter(1,2);
ACGiter(3,3)=summary_tables.iter(1,3);

filename='104Vector14instances.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');

%%%%%% Problem 4
% Create global hyperparams
N = 1000;
seed = 91;
dimM = 20;
dimN = 1000;
M=10^1;
m=10^1;
global_tol = 1e-4;
time_limit = 10800;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

[oracle, hparams] = test_fn_lin_cone_constr_01(N, M, m, seed, dimM, dimN);

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
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL','RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);


FVAL(4,1)=summary_tables.fval(1,1);
FVAL(4,2)=summary_tables.fval(1,2);
FVAL(4,3)=summary_tables.fval(1,3);

RUN(4,1)=summary_tables.runtime(1,1);
RUN(4,2)=summary_tables.runtime(1,2);
RUN(4,3)=summary_tables.runtime(1,3);


ACGiter(4,1)=summary_tables.iter(1,1);
ACGiter(4,2)=summary_tables.iter(1,2);
ACGiter(4,3)=summary_tables.iter(1,3);


filename='104Vector14instances.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');


%%%%% Problem 5
% Create global hyperparams
N = 1000;
seed = 556;
dimM = 20;
dimN = 1000;
M=10^2;
m=10^1;
global_tol = 1e-4;
time_limit = 10800;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

[oracle, hparams] = test_fn_lin_cone_constr_01(N, M, m, seed, dimM, dimN);

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
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL','RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);


FVAL(5,1)=summary_tables.fval(1,1);
FVAL(5,2)=summary_tables.fval(1,2);
FVAL(5,3)=summary_tables.fval(1,3);

RUN(5,1)=summary_tables.runtime(1,1);
RUN(5,2)=summary_tables.runtime(1,2);
RUN(5,3)=summary_tables.runtime(1,3);


ACGiter(5,1)=summary_tables.iter(1,1);
ACGiter(5,2)=summary_tables.iter(1,2);
ACGiter(5,3)=summary_tables.iter(1,3);


filename='104Vector14instances.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');


%%%%%%% Problem 6
% Create global hyperparams
N = 1000;
seed = 224;
dimM = 20;
dimN = 1000;
M=10^3;
m=10^1;
global_tol = 1e-4;
time_limit = 10800;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

[oracle, hparams] = test_fn_lin_cone_constr_01(N, M, m, seed, dimM, dimN);

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
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL','RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);


FVAL(6,1)=summary_tables.fval(1,1);
FVAL(6,2)=summary_tables.fval(1,2);
FVAL(6,3)=summary_tables.fval(1,3);

RUN(6,1)=summary_tables.runtime(1,1);
RUN(6,2)=summary_tables.runtime(1,2);
RUN(6,3)=summary_tables.runtime(1,3);


ACGiter(6,1)=summary_tables.iter(1,1);
ACGiter(6,2)=summary_tables.iter(1,2);
ACGiter(6,3)=summary_tables.iter(1,3);


filename='104Vector14instances.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');

%%%%%%% Problem 7
% Create global hyperparams
N = 1000;
seed = 392;
dimM = 20;
dimN = 1000;
M=10^4;
m=10^1;
global_tol = 1e-4;
time_limit = 10800;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

[oracle, hparams] = test_fn_lin_cone_constr_01(N, M, m, seed, dimM, dimN);

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
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL','RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);


FVAL(7,1)=summary_tables.fval(1,1);
FVAL(7,2)=summary_tables.fval(1,2);
FVAL(7,3)=summary_tables.fval(1,3);

RUN(7,1)=summary_tables.runtime(1,1);
RUN(7,2)=summary_tables.runtime(1,2);
RUN(7,3)=summary_tables.runtime(1,3);

ACGiter(7,1)=summary_tables.iter(1,1);
ACGiter(7,2)=summary_tables.iter(1,2);
ACGiter(7,3)=summary_tables.iter(1,3);


filename='104Vector14instances.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');


%%%% Problem 8
% Create global hyperparams
N = 1000;
seed = 575;
dimM = 20;
dimN = 1000;
M=10^3;
m=10^2;
global_tol = 1e-4;
time_limit = 10800;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

[oracle, hparams] = test_fn_lin_cone_constr_01(N, M, m, seed, dimM, dimN);

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
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL','RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);


FVAL(8,1)=summary_tables.fval(1,1);
FVAL(8,2)=summary_tables.fval(1,2);
FVAL(8,3)=summary_tables.fval(1,3);

RUN(8,1)=summary_tables.runtime(1,1);
RUN(8,2)=summary_tables.runtime(1,2);
RUN(8,3)=summary_tables.runtime(1,3);

ACGiter(8,1)=summary_tables.iter(1,1);
ACGiter(8,2)=summary_tables.iter(1,2);
ACGiter(8,3)=summary_tables.iter(1,3);

filename='104Vector14instances.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');

%%%% Problem 9
% Create global hyperparams
N = 1000;
seed = 614;
dimM = 20;
dimN = 1000;
M=10^4;
m=10^2;
global_tol = 1e-4;
time_limit = 10800;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

[oracle, hparams] = test_fn_lin_cone_constr_01(N, M, m, seed, dimM, dimN);

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
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL','RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);


FVAL(9,1)=summary_tables.fval(1,1);
FVAL(9,2)=summary_tables.fval(1,2);
FVAL(9,3)=summary_tables.fval(1,3);

RUN(9,1)=summary_tables.runtime(1,1);
RUN(9,2)=summary_tables.runtime(1,2);
RUN(9,3)=summary_tables.runtime(1,3);

ACGiter(9,1)=summary_tables.iter(1,1);
ACGiter(9,2)=summary_tables.iter(1,2);
ACGiter(9,3)=summary_tables.iter(1,3);


filename='104Vector14instances.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');

%%%%% Problem 10
% Create global hyperparams
N = 1000;
seed = 442;
dimM = 20;
dimN = 1000;
M=10^5;
m=10^2;
global_tol = 1e-4;
time_limit = 10800;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

[oracle, hparams] = test_fn_lin_cone_constr_01(N, M, m, seed, dimM, dimN);

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
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL','RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);


FVAL(10,1)=summary_tables.fval(1,1);
FVAL(10,2)=summary_tables.fval(1,2);
FVAL(10,3)=summary_tables.fval(1,3);

RUN(10,1)=summary_tables.runtime(1,1);
RUN(10,2)=summary_tables.runtime(1,2);
RUN(10,3)=summary_tables.runtime(1,3);

ACGiter(10,1)=summary_tables.iter(1,1);
ACGiter(10,2)=summary_tables.iter(1,2);
ACGiter(10,3)=summary_tables.iter(1,3);

filename='104Vector14instances.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');

%%%%%%% Problem 11
% Create global hyperparams
N = 1000;
seed = 135;
dimM = 20;
dimN = 1000;
M=10^3;
m=10^3;
global_tol = 1e-4;
time_limit = 10800;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

[oracle, hparams] = test_fn_lin_cone_constr_01(N, M, m, seed, dimM, dimN);

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
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL','RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);


FVAL(11,1)=summary_tables.fval(1,1);
FVAL(11,2)=summary_tables.fval(1,2);
FVAL(11,3)=summary_tables.fval(1,3);

RUN(11,1)=summary_tables.runtime(1,1);
RUN(11,2)=summary_tables.runtime(1,2);
RUN(11,3)=summary_tables.runtime(1,3);

ACGiter(11,1)=summary_tables.iter(1,1);
ACGiter(11,2)=summary_tables.iter(1,2);
ACGiter(11,3)=summary_tables.iter(1,3);


filename='104Vector14instances.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');

%%%% Problem 12
% Create global hyperparams
N = 1000;
seed = 177;
dimM = 20;
dimN = 1000;
M=10^4;
m=10^3;
global_tol = 1e-4;
time_limit = 10800;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

[oracle, hparams] = test_fn_lin_cone_constr_01(N, M, m, seed, dimM, dimN);

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
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL','RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);


FVAL(12,1)=summary_tables.fval(1,1);
FVAL(12,2)=summary_tables.fval(1,2);
FVAL(12,3)=summary_tables.fval(1,3);

RUN(12,1)=summary_tables.runtime(1,1);
RUN(12,2)=summary_tables.runtime(1,2);
RUN(12,3)=summary_tables.runtime(1,3);


ACGiter(12,1)=summary_tables.iter(1,1);
ACGiter(12,2)=summary_tables.iter(1,2);
ACGiter(12,3)=summary_tables.iter(1,3);


filename='104Vector14instances.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');

%%% Problem 13

% Create global hyperparams
N = 1000;
seed = 291;
dimM = 20;
dimN = 1000;
M=10^5;
m=10^3;
global_tol = 1e-4;
time_limit = 10800;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

[oracle, hparams] = test_fn_lin_cone_constr_01(N, M, m, seed, dimM, dimN);

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
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL','RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);


FVAL(13,1)=summary_tables.fval(1,1);
FVAL(13,2)=summary_tables.fval(1,2);
FVAL(13,3)=summary_tables.fval(1,3);

RUN(13,1)=summary_tables.runtime(1,1);
RUN(13,2)=summary_tables.runtime(1,2);
RUN(13,3)=summary_tables.runtime(1,3);


ACGiter(13,1)=summary_tables.iter(1,1);
ACGiter(13,2)=summary_tables.iter(1,2);
ACGiter(13,3)=summary_tables.iter(1,3);

filename='104Vector14instances.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');


%%%%%% Problem 14
% Create global hyperparams
N = 1000;
seed = 714;
dimM = 20;
dimN = 1000;
M=10^6;
m=10^3;
global_tol = 1e-4;
time_limit = 10800;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

[oracle, hparams] = test_fn_lin_cone_constr_01(N, M, m, seed, dimM, dimN);

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
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL','RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(ncvx_lc_qp, framework_arr, solver_arr, hparam_arr, name_arr);
disp(summary_tables.all);


FVAL(14,1)=summary_tables.fval(1,1);
FVAL(14,2)=summary_tables.fval(1,2);
FVAL(14,3)=summary_tables.fval(1,3);

RUN(14,1)=summary_tables.runtime(1,1);
RUN(14,2)=summary_tables.runtime(1,2);
RUN(14,3)=summary_tables.runtime(1,3);


ACGiter(14,1)=summary_tables.iter(1,1);
ACGiter(14,2)=summary_tables.iter(1,2);
ACGiter(14,3)=summary_tables.iter(1,3);

filename='104Vector14instances.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A21');
writetable(ACGiter,filename,'Sheet',1,'Range','A41');