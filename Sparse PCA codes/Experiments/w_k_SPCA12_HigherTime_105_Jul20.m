%% SPDX-License-Identifier: MIT
% Copyright © 2021 Weiwei "William" Kong

% Solve a sparse PCA problem.


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

%%%%%Problem 1
Problem=1
% Create global hyperparams
b = 0.008;
nu = 100;
p = 100;
n = 100;
k = 5;
s=5;
seed = 648;
global_tol = 1e-5;
time_limit = 3600;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_spca_01(b, nu, p, n, s, k, seed);

% Problem dependent hparams
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Create the Model object and specify the solver.
spca = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
spca.M = hparams.M;
spca.m = hparams.m;
spca.x0 = hparams.x0;
spca.K_constr = hparams.K_constr;

% Add linear constraints.
spca.constr_fn = @(x) hparams.constr_fn(x);
spca.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
spca.opt_type = 'relative';
spca.feas_type = 'relative';
spca.opt_tol = global_tol;
spca.feas_tol = global_tol;
spca.time_limit = time_limit;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL', 'RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(spca, framework_arr, solver_arr, hparam_arr, name_arr);

disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);
FVAL(Problem,3)=summary_tables.fval(1,3);

RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);
RUN(Problem,3)=summary_tables.runtime(1,3);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);
ACGiter(Problem,3)=summary_tables.iter(1,3);

filename='105SPCA_HigherTime_12inst.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A31');
writetable(ACGiter,filename,'Sheet',1,'Range','A61');

%%%%%Problem 2
Problem=2
% Create global hyperparams
b = 0.008;
nu = 100;
p = 100;
n = 100;
k = 10;
s=5;
seed = 648;
global_tol = 1e-5;
time_limit = 3600;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_spca_01(b, nu, p, n, s, k, seed);

% Problem dependent hparams
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Create the Model object and specify the solver.
spca = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
spca.M = hparams.M;
spca.m = hparams.m;
spca.x0 = hparams.x0;
spca.K_constr = hparams.K_constr;

% Add linear constraints.
spca.constr_fn = @(x) hparams.constr_fn(x);
spca.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
spca.opt_type = 'relative';
spca.feas_type = 'relative';
spca.opt_tol = global_tol;
spca.feas_tol = global_tol;
spca.time_limit = time_limit;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL', 'RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(spca, framework_arr, solver_arr, hparam_arr, name_arr);

disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);
FVAL(Problem,3)=summary_tables.fval(1,3);

RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);
RUN(Problem,3)=summary_tables.runtime(1,3);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);
ACGiter(Problem,3)=summary_tables.iter(1,3);

filename='105SPCA_HigherTime_12inst.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A31');
writetable(ACGiter,filename,'Sheet',1,'Range','A61');

%%%%%Problem 3
Problem=3
% Create global hyperparams
b = 0.008;
nu = 100;
p = 100;
n = 100;
k = 20;
s=5;
seed = 648;
global_tol = 1e-5;
time_limit = 3600;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_spca_01(b, nu, p, n, s, k, seed);

% Problem dependent hparams
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Create the Model object and specify the solver.
spca = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
spca.M = hparams.M;
spca.m = hparams.m;
spca.x0 = hparams.x0;
spca.K_constr = hparams.K_constr;

% Add linear constraints.
spca.constr_fn = @(x) hparams.constr_fn(x);
spca.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
spca.opt_type = 'relative';
spca.feas_type = 'relative';
spca.opt_tol = global_tol;
spca.feas_tol = global_tol;
spca.time_limit = time_limit;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL', 'RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(spca, framework_arr, solver_arr, hparam_arr, name_arr);

disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);
FVAL(Problem,3)=summary_tables.fval(1,3);

RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);
RUN(Problem,3)=summary_tables.runtime(1,3);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);
ACGiter(Problem,3)=summary_tables.iter(1,3);

filename='105SPCA_HigherTime_12inst.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A31');
writetable(ACGiter,filename,'Sheet',1,'Range','A61');


%%%%%Problem 4
Problem=4
% Create global hyperparams
b = 0.005;
nu = 100;
p = 100;
n = 100;
k = 5;
s=5;
seed = 442;
global_tol = 1e-5;
time_limit = 3600;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_spca_01(b, nu, p, n, s, k, seed);

% Problem dependent hparams
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Create the Model object and specify the solver.
spca = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
spca.M = hparams.M;
spca.m = hparams.m;
spca.x0 = hparams.x0;
spca.K_constr = hparams.K_constr;

% Add linear constraints.
spca.constr_fn = @(x) hparams.constr_fn(x);
spca.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
spca.opt_type = 'relative';
spca.feas_type = 'relative';
spca.opt_tol = global_tol;
spca.feas_tol = global_tol;
spca.time_limit = time_limit;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL', 'RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(spca, framework_arr, solver_arr, hparam_arr, name_arr);

disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);
FVAL(Problem,3)=summary_tables.fval(1,3);

RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);
RUN(Problem,3)=summary_tables.runtime(1,3);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);
ACGiter(Problem,3)=summary_tables.iter(1,3);

filename='105SPCA_HigherTime_12inst.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A31');
writetable(ACGiter,filename,'Sheet',1,'Range','A61');


%%%%%Problem 5
Problem=5
% Create global hyperparams
b = 0.005;
nu = 100;
p = 100;
n = 100;
k = 10;
s=5;
seed = 442;
global_tol = 1e-5;
time_limit = 3600;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_spca_01(b, nu, p, n, s, k, seed);

% Problem dependent hparams
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Create the Model object and specify the solver.
spca = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
spca.M = hparams.M;
spca.m = hparams.m;
spca.x0 = hparams.x0;
spca.K_constr = hparams.K_constr;

% Add linear constraints.
spca.constr_fn = @(x) hparams.constr_fn(x);
spca.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
spca.opt_type = 'relative';
spca.feas_type = 'relative';
spca.opt_tol = global_tol;
spca.feas_tol = global_tol;
spca.time_limit = time_limit;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL', 'RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(spca, framework_arr, solver_arr, hparam_arr, name_arr);

disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);
FVAL(Problem,3)=summary_tables.fval(1,3);

RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);
RUN(Problem,3)=summary_tables.runtime(1,3);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);
ACGiter(Problem,3)=summary_tables.iter(1,3);

filename='105SPCA_HigherTime_12inst.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A31');
writetable(ACGiter,filename,'Sheet',1,'Range','A61');

%%%%%Problem 6
Problem=6
% Create global hyperparams
b = 0.005;
nu = 100;
p = 100;
n = 100;
k = 20;
s=5;
seed = 442;
global_tol = 1e-5;
time_limit = 3600;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_spca_01(b, nu, p, n, s, k, seed);

% Problem dependent hparams
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Create the Model object and specify the solver.
spca = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
spca.M = hparams.M;
spca.m = hparams.m;
spca.x0 = hparams.x0;
spca.K_constr = hparams.K_constr;

% Add linear constraints.
spca.constr_fn = @(x) hparams.constr_fn(x);
spca.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
spca.opt_type = 'relative';
spca.feas_type = 'relative';
spca.opt_tol = global_tol;
spca.feas_tol = global_tol;
spca.time_limit = time_limit;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL', 'RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(spca, framework_arr, solver_arr, hparam_arr, name_arr);

disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);
FVAL(Problem,3)=summary_tables.fval(1,3);

RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);
RUN(Problem,3)=summary_tables.runtime(1,3);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);
ACGiter(Problem,3)=summary_tables.iter(1,3);

filename='105SPCA_HigherTime_12inst.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A31');
writetable(ACGiter,filename,'Sheet',1,'Range','A61');

%%%%%Problem 7
Problem=7
% Create global hyperparams
b = 0.004;
nu = 100;
p = 100;
n = 100;
k = 5;
s=5;
seed = 512;
global_tol = 1e-5;
time_limit = 3600;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_spca_01(b, nu, p, n, s, k, seed);

% Problem dependent hparams
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Create the Model object and specify the solver.
spca = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
spca.M = hparams.M;
spca.m = hparams.m;
spca.x0 = hparams.x0;
spca.K_constr = hparams.K_constr;

% Add linear constraints.
spca.constr_fn = @(x) hparams.constr_fn(x);
spca.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
spca.opt_type = 'relative';
spca.feas_type = 'relative';
spca.opt_tol = global_tol;
spca.feas_tol = global_tol;
spca.time_limit = time_limit;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL', 'RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(spca, framework_arr, solver_arr, hparam_arr, name_arr);

disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);
FVAL(Problem,3)=summary_tables.fval(1,3);

RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);
RUN(Problem,3)=summary_tables.runtime(1,3);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);
ACGiter(Problem,3)=summary_tables.iter(1,3);

filename='105SPCA_HigherTime_12inst.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A31');
writetable(ACGiter,filename,'Sheet',1,'Range','A61');


%%%%%Problem 8
Problem=8
% Create global hyperparams
b = 0.004;
nu = 100;
p = 100;
n = 100;
k = 10;
s=5;
seed = 512;
global_tol = 1e-5;
time_limit = 3600;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_spca_01(b, nu, p, n, s, k, seed);

% Problem dependent hparams
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Create the Model object and specify the solver.
spca = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
spca.M = hparams.M;
spca.m = hparams.m;
spca.x0 = hparams.x0;
spca.K_constr = hparams.K_constr;

% Add linear constraints.
spca.constr_fn = @(x) hparams.constr_fn(x);
spca.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
spca.opt_type = 'relative';
spca.feas_type = 'relative';
spca.opt_tol = global_tol;
spca.feas_tol = global_tol;
spca.time_limit = time_limit;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL', 'RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(spca, framework_arr, solver_arr, hparam_arr, name_arr);

disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);
FVAL(Problem,3)=summary_tables.fval(1,3);

RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);
RUN(Problem,3)=summary_tables.runtime(1,3);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);
ACGiter(Problem,3)=summary_tables.iter(1,3);

filename='105SPCA_HigherTime_12inst.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A31');
writetable(ACGiter,filename,'Sheet',1,'Range','A61');

%%%%%Problem 9
Problem=9
% Create global hyperparams
b = 0.004;
nu = 100;
p = 100;
n = 100;
k = 20;
s=5;
seed = 512;
global_tol = 1e-5;
time_limit = 3600;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_spca_01(b, nu, p, n, s, k, seed);

% Problem dependent hparams
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Create the Model object and specify the solver.
spca = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
spca.M = hparams.M;
spca.m = hparams.m;
spca.x0 = hparams.x0;
spca.K_constr = hparams.K_constr;

% Add linear constraints.
spca.constr_fn = @(x) hparams.constr_fn(x);
spca.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
spca.opt_type = 'relative';
spca.feas_type = 'relative';
spca.opt_tol = global_tol;
spca.feas_tol = global_tol;
spca.time_limit = time_limit;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL', 'RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(spca, framework_arr, solver_arr, hparam_arr, name_arr);

disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);
FVAL(Problem,3)=summary_tables.fval(1,3);

RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);
RUN(Problem,3)=summary_tables.runtime(1,3);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);
ACGiter(Problem,3)=summary_tables.iter(1,3);

filename='105SPCA_HigherTime_12inst.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A31');
writetable(ACGiter,filename,'Sheet',1,'Range','A61');

%%%%%Problem 10
Problem=10
% Create global hyperparams
b = 0.001;
nu = 100;
p = 100;
n = 100;
k = 5;
s=5;
seed = 301;
global_tol = 1e-5;
time_limit = 3600;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_spca_01(b, nu, p, n, s, k, seed);

% Problem dependent hparams
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Create the Model object and specify the solver.
spca = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
spca.M = hparams.M;
spca.m = hparams.m;
spca.x0 = hparams.x0;
spca.K_constr = hparams.K_constr;

% Add linear constraints.
spca.constr_fn = @(x) hparams.constr_fn(x);
spca.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
spca.opt_type = 'relative';
spca.feas_type = 'relative';
spca.opt_tol = global_tol;
spca.feas_tol = global_tol;
spca.time_limit = time_limit;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL', 'RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(spca, framework_arr, solver_arr, hparam_arr, name_arr);

disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);
FVAL(Problem,3)=summary_tables.fval(1,3);

RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);
RUN(Problem,3)=summary_tables.runtime(1,3);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);
ACGiter(Problem,3)=summary_tables.iter(1,3);

filename='105SPCA_HigherTime_12inst.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A31');
writetable(ACGiter,filename,'Sheet',1,'Range','A61');

%%%%%Problem 11
Problem=11
% Create global hyperparams
b = 0.001;
nu = 100;
p = 100;
n = 100;
k = 10;
s=5;
seed = 301;
global_tol = 1e-5;
time_limit = 3600;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_spca_01(b, nu, p, n, s, k, seed);

% Problem dependent hparams
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Create the Model object and specify the solver.
spca = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
spca.M = hparams.M;
spca.m = hparams.m;
spca.x0 = hparams.x0;
spca.K_constr = hparams.K_constr;

% Add linear constraints.
spca.constr_fn = @(x) hparams.constr_fn(x);
spca.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
spca.opt_type = 'relative';
spca.feas_type = 'relative';
spca.opt_tol = global_tol;
spca.feas_tol = global_tol;
spca.time_limit = time_limit;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL', 'RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(spca, framework_arr, solver_arr, hparam_arr, name_arr);

disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);
FVAL(Problem,3)=summary_tables.fval(1,3);

RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);
RUN(Problem,3)=summary_tables.runtime(1,3);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);
ACGiter(Problem,3)=summary_tables.iter(1,3);

filename='105SPCA_HigherTime_12inst.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A31');
writetable(ACGiter,filename,'Sheet',1,'Range','A61');

%%%%%Problem 12
Problem=12
% Create global hyperparams
b = 0.001;
nu = 100;
p = 100;
n = 100;
k = 20;
s=5;
seed = 301;
global_tol = 1e-5;
time_limit = 3600;

% -------------------------------------------------------------------------
%% Table 1
% -------------------------------------------------------------------------
disp('========')
disp('TABLE 1');
disp('========')

% Use a problem instance generator to create the oracle and
% hyperparameters.

[oracle, hparams] = test_fn_spca_01(b, nu, p, n, s, k, seed);

% Problem dependent hparams
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
% Note that we are using the fact that |X|_F <= 1 over the spectraplex.
ialm_hparam.B_vec = hparams.K_constr_vec;

% Create the Model object and specify the solver.
spca = ConstrCompModel(oracle);

% Set the curvatures, the starting point x0, and special functions.
spca.M = hparams.M;
spca.m = hparams.m;
spca.x0 = hparams.x0;
spca.K_constr = hparams.K_constr;

% Add linear constraints.
spca.constr_fn = @(x) hparams.constr_fn(x);
spca.grad_constr_fn = hparams.grad_constr_fn;

% Set up the termination criterion.
spca.opt_type = 'relative';
spca.feas_type = 'relative';
spca.opt_tol = global_tol;
spca.feas_tol = global_tol;
spca.time_limit = time_limit;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam};
name_arr = {'iALM', 'IPL', 'RQP'};
framework_arr = {@iALM3, @IAIPAL3, @penalty3};
solver_arr = {@ECG, @ECG, @AIPP};
[summary_tables, comp_models] = run_CCM_benchmark(spca, framework_arr, solver_arr, hparam_arr, name_arr);

disp(summary_tables.all);

FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);
FVAL(Problem,3)=summary_tables.fval(1,3);

RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);
RUN(Problem,3)=summary_tables.runtime(1,3);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);
ACGiter(Problem,3)=summary_tables.iter(1,3);

filename='105SPCA_HigherTime_12inst.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A31');
writetable(ACGiter,filename,'Sheet',1,'Range','A61');