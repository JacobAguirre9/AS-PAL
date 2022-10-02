%% SPDX-License-Identifier: MIT
% Copyright © 2021 Weiwei "William" Kong

% Solve a multivariate nonconvex linearly constrained quadratic programming problem constrained to a box.


format long
%%%Problem=1
% Initialize
Problem=1
N = 1000;
r=5;
seed=976;
dimM=20;
dimN=100;
M=10^1;
m=10^0;
global_tol = 1e-5;



[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();

%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;
%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;

%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {iapial_hparam};
name_arr = {'IPL'};
framework_arr = {@IAIPAL3};
solver_arr = {@ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
RUN(Problem,1)=summary_tables.runtime(1,1);
ACGiter(Problem,1)=summary_tables.iter(1,1);


filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');


%%%Problem=2
% Initialize
Problem=2
N = 1000;
r=10;
seed=976;
dimM=20;
dimN=100;
M=10^1;
m=10^0;
global_tol = 1e-5;



[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();

%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';
rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;
%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;





%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {iapial_hparam};
name_arr = {'IPL'};
framework_arr = {@IAIPAL3};
solver_arr = {@ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
RUN(Problem,1)=summary_tables.runtime(1,1);
ACGiter(Problem,1)=summary_tables.iter(1,1);


filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');

%%%Problem=3
% Initialize
Problem=3
N = 1000;
r=20;
seed=976;
dimM=20;
dimN=100;
M=10^1;
m=10^0;
global_tol = 1e-5;

[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();

%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;

%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;

%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {iapial_hparam};
name_arr = {'IPL'};
framework_arr = {@IAIPAL3};
solver_arr = {@ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
RUN(Problem,1)=summary_tables.runtime(1,1);
ACGiter(Problem,1)=summary_tables.iter(1,1);

filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');


%%%Problem=4
% Initialize
Problem=4
N = 1000;
r=1;
seed=167;
dimM=20;
dimN=100;
M=10^2;
m=10^1;
global_tol = 1e-5;



[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();

%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;

%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;

%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {iapial_hparam};
name_arr = {'IPL'};
framework_arr = {@IAIPAL3};
solver_arr = {@ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
RUN(Problem,1)=summary_tables.runtime(1,1);
ACGiter(Problem,1)=summary_tables.iter(1,1);


filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');



%%%Problem=5
% Initialize
Problem=5
N = 1000;
r=2;
seed=167;
dimM=20;
dimN=100;
M=10^2;
m=10^1;
global_tol = 1e-5;

[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();

%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;

%%%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;

%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {iapial_hparam};
name_arr = {'IPL'};
framework_arr = {@IAIPAL3};
solver_arr = {@ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
RUN(Problem,1)=summary_tables.runtime(1,1);
ACGiter(Problem,1)=summary_tables.iter(1,1);


filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');

%%%Problem=6
% Initialize
Problem=6
N = 1000;
r=5;
seed=167;
dimM=20;
dimN=100;
M=10^2;
m=10^1;
global_tol = 1e-5;

[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();

%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;

%%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;


%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {iapial_hparam};
name_arr = {'IPL'};
framework_arr = {@IAIPAL3};
solver_arr = {@ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
RUN(Problem,1)=summary_tables.runtime(1,1);
ACGiter(Problem,1)=summary_tables.iter(1,1);


filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');



%%%Problem=7
% Initialize
Problem=7
N = 1000;
r=1;
seed=718;
dimM=20;
dimN=100;
M=10^3;
m=10^1;
global_tol = 1e-5;

[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();


%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;
%%%%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;

%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {iapial_hparam};
name_arr = {'IPL'};
framework_arr = {@IAIPAL3};
solver_arr = {@ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
RUN(Problem,1)=summary_tables.runtime(1,1);
ACGiter(Problem,1)=summary_tables.iter(1,1);


filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');


%%%Problem=8
% Initialize
Problem=8
N = 1000;
r=2;
seed=718;
dimM=20;
dimN=100;
M=10^3;
m=10^1;
global_tol = 1e-5;

[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();


%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;
%%%%%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;

%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {iapial_hparam};
name_arr = {'IPL'};
framework_arr = {@IAIPAL3};
solver_arr = {@ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
RUN(Problem,1)=summary_tables.runtime(1,1);
ACGiter(Problem,1)=summary_tables.iter(1,1);

filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');

%%%Problem=9
% Initialize
Problem=9
N = 1000;
r=5;
seed=718;
dimM=20;
dimN=100;
M=10^3;
m=10^1;
global_tol = 1e-5;

[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();


%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;
%%%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;

%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {iapial_hparam};
name_arr = {'IPL'};
framework_arr = {@IAIPAL3};
solver_arr = {@ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
RUN(Problem,1)=summary_tables.runtime(1,1);
ACGiter(Problem,1)=summary_tables.iter(1,1);

filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');


%%%Problem=10
% Initialize
Problem=10
N = 1000;
r=1;
seed=771;
dimM=20;
dimN=100;
M=10^3;
m=10^2;
global_tol = 1e-5;

[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();


%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;
%%%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;
%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {iapial_hparam};
name_arr = {'IPL'};
framework_arr = {@IAIPAL3};
solver_arr = {@ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
RUN(Problem,1)=summary_tables.runtime(1,1);
ACGiter(Problem,1)=summary_tables.iter(1,1);

filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');


%%%Problem=11
% Initialize
Problem=11
N = 1000;
r=2;
seed=771;
dimM=20;
dimN=100;
M=10^3;
m=10^2;
global_tol = 1e-5;

[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();


%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;
%%%%%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;
%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {iapial_hparam};
name_arr = {'IPL'};
framework_arr = {@IAIPAL3};
solver_arr = {@ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
RUN(Problem,1)=summary_tables.runtime(1,1);
ACGiter(Problem,1)=summary_tables.iter(1,1);

filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');


%%%Problem=12
% Initialize
Problem=12
N = 1000;
r=5;
seed=771;
dimM=20;
dimN=100;
M=10^3;
m=10^2;
global_tol = 1e-5;

[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();


%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;
%%%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;

%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {iapial_hparam};
name_arr = {'IPL'};
framework_arr = {@IAIPAL3};
solver_arr = {@ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
RUN(Problem,1)=summary_tables.runtime(1,1);
ACGiter(Problem,1)=summary_tables.iter(1,1);

filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');

%%%Problem=13
% Initialize
Problem=13
N = 1000;
r=1;
seed=418;
dimM=20;
dimN=100;
M=10^4;
m=10^2;
global_tol = 1e-5;

[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();


%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;
%%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;

%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam, spa1_hparam, spa2_hparam}; 
name_arr = {'iALM', 'IPL', 'RQP', 'SPA1', 'SPA2'};
framework_arr = {@iALM3, @IAIPAL3,  @penalty3, @sProxALM, @sProxALM};
solver_arr = {@ECG, @ECG, @AIPP, @ECG, @ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);
FVAL(Problem,3)=summary_tables.fval(1,3);
FVAL(Problem,4)=summary_tables.fval(1,4);
FVAL(Problem,5)=summary_tables.fval(1,5);

RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);
RUN(Problem,3)=summary_tables.runtime(1,3);
RUN(Problem,4)=summary_tables.runtime(1,4);
RUN(Problem,5)=summary_tables.runtime(1,5);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);
ACGiter(Problem,3)=summary_tables.iter(1,3);
ACGiter(Problem,4)=summary_tables.iter(1,4);
ACGiter(Problem,5)=summary_tables.iter(1,5);

filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');


%%%Problem=14
% Initialize
Problem=14
N = 1000;
r=2;
seed=418;
dimM=20;
dimN=100;
M=10^4;
m=10^2;
global_tol = 1e-5;

[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();


%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;
%%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;

%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam, spa1_hparam, spa2_hparam}; 
name_arr = {'iALM', 'IPL', 'RQP', 'SPA1', 'SPA2'};
framework_arr = {@iALM3, @IAIPAL3,  @penalty3, @sProxALM, @sProxALM};
solver_arr = {@ECG, @ECG, @AIPP, @ECG, @ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);
FVAL(Problem,3)=summary_tables.fval(1,3);
FVAL(Problem,4)=summary_tables.fval(1,4);
FVAL(Problem,5)=summary_tables.fval(1,5);

RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);
RUN(Problem,3)=summary_tables.runtime(1,3);
RUN(Problem,4)=summary_tables.runtime(1,4);
RUN(Problem,5)=summary_tables.runtime(1,5);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);
ACGiter(Problem,3)=summary_tables.iter(1,3);
ACGiter(Problem,4)=summary_tables.iter(1,4);
ACGiter(Problem,5)=summary_tables.iter(1,5);

filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');



%%%Problem=15
% Initialize
Problem=15
N = 1000;
r=5;
seed=418;
dimM=20;
dimN=100;
M=10^4;
m=10^2;
global_tol = 1e-5;

[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();


%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;
%%%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;

%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {iapial_hparam};
name_arr = {'IPL'};
framework_arr = {@IAIPAL3};
solver_arr = {@ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
RUN(Problem,1)=summary_tables.runtime(1,1);
ACGiter(Problem,1)=summary_tables.iter(1,1);


filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');


%%%Problem=16
% Initialize
Problem=16
N = 1000;
r=5;
seed=323;
dimM=20;
dimN=100;
M=10^3;
m=10^3;
global_tol = 1e-5;

[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();


%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;
%%%%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;

%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {iapial_hparam};
name_arr = {'IPL'};
framework_arr = {@IAIPAL3};
solver_arr = {@ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
RUN(Problem,1)=summary_tables.runtime(1,1);
ACGiter(Problem,1)=summary_tables.iter(1,1);

filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');


%%%Problem=17
% Initialize
Problem=17
N = 1000;
r=10;
seed=323;
dimM=20;
dimN=100;
M=10^3;
m=10^3;
global_tol = 1e-5;

[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();


%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;
%%%%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;

%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {iapial_hparam};
name_arr = {'IPL'};
framework_arr = {@IAIPAL3};
solver_arr = {@ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
RUN(Problem,1)=summary_tables.runtime(1,1);
ACGiter(Problem,1)=summary_tables.iter(1,1);

filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');


%%%Problem=18
% Initialize
Problem=18
N = 1000;
r=1;
seed=887;
dimM=20;
dimN=100;
M=10^4;
m=10^3;
global_tol = 1e-5;

[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();

%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;
%%%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;

%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam, spa1_hparam, spa2_hparam}; 
name_arr = {'iALM', 'IPL', 'RQP', 'SPA1', 'SPA2'};
framework_arr = {@iALM3, @IAIPAL3,  @penalty3, @sProxALM, @sProxALM};
solver_arr = {@ECG, @ECG, @AIPP, @ECG, @ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);
FVAL(Problem,3)=summary_tables.fval(1,3);
FVAL(Problem,4)=summary_tables.fval(1,4);
FVAL(Problem,5)=summary_tables.fval(1,5);

RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);
RUN(Problem,3)=summary_tables.runtime(1,3);
RUN(Problem,4)=summary_tables.runtime(1,4);
RUN(Problem,5)=summary_tables.runtime(1,5);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);
ACGiter(Problem,3)=summary_tables.iter(1,3);
ACGiter(Problem,4)=summary_tables.iter(1,4);
ACGiter(Problem,5)=summary_tables.iter(1,5);

filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');



%%%Problem=19
% Initialize
Problem=19
N = 1000;
r=2;
seed=887;
dimM=20;
dimN=100;
M=10^4;
m=10^3;
global_tol = 1e-5;

[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();


%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;
%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;
%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam, spa1_hparam, spa2_hparam}; 
name_arr = {'iALM', 'IPL', 'RQP', 'SPA1', 'SPA2'};
framework_arr = {@iALM3, @IAIPAL3,  @penalty3, @sProxALM, @sProxALM};
solver_arr = {@ECG, @ECG, @AIPP, @ECG, @ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);
FVAL(Problem,3)=summary_tables.fval(1,3);
FVAL(Problem,4)=summary_tables.fval(1,4);
FVAL(Problem,5)=summary_tables.fval(1,5);

RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);
RUN(Problem,3)=summary_tables.runtime(1,3);
RUN(Problem,4)=summary_tables.runtime(1,4);
RUN(Problem,5)=summary_tables.runtime(1,5);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);
ACGiter(Problem,3)=summary_tables.iter(1,3);
ACGiter(Problem,4)=summary_tables.iter(1,4);
ACGiter(Problem,5)=summary_tables.iter(1,5);

filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');


%%%Problem=20
% Initialize
Problem=20
N = 1000;
r=5;
seed=887;
dimM=20;
dimN=100;
M=10^4;
m=10^3;
global_tol = 1e-5;

[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();


%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;
%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;

%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam, spa1_hparam, spa2_hparam}; 
name_arr = {'iALM', 'IPL', 'RQP', 'SPA1', 'SPA2'};
framework_arr = {@iALM3, @IAIPAL3,  @penalty3, @sProxALM, @sProxALM};
solver_arr = {@ECG, @ECG, @AIPP, @ECG, @ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);
FVAL(Problem,3)=summary_tables.fval(1,3);
FVAL(Problem,4)=summary_tables.fval(1,4);
FVAL(Problem,5)=summary_tables.fval(1,5);

RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);
RUN(Problem,3)=summary_tables.runtime(1,3);
RUN(Problem,4)=summary_tables.runtime(1,4);
RUN(Problem,5)=summary_tables.runtime(1,5);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);
ACGiter(Problem,3)=summary_tables.iter(1,3);
ACGiter(Problem,4)=summary_tables.iter(1,4);
ACGiter(Problem,5)=summary_tables.iter(1,5);

filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');


%%%Problem=21
% Initialize
Problem=21
N = 1000;
r=1;
seed=218;
dimM=20;
dimN=100;
M=10^5;
m=10^3;
global_tol = 1e-5;

[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();


%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;
%%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;

%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam, spa1_hparam, spa2_hparam}; 
name_arr = {'iALM', 'IPL', 'RQP', 'SPA1', 'SPA2'};
framework_arr = {@iALM3, @IAIPAL3,  @penalty3, @sProxALM, @sProxALM};
solver_arr = {@ECG, @ECG, @AIPP, @ECG, @ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);
FVAL(Problem,3)=summary_tables.fval(1,3);
FVAL(Problem,4)=summary_tables.fval(1,4);
FVAL(Problem,5)=summary_tables.fval(1,5);

RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);
RUN(Problem,3)=summary_tables.runtime(1,3);
RUN(Problem,4)=summary_tables.runtime(1,4);
RUN(Problem,5)=summary_tables.runtime(1,5);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);
ACGiter(Problem,3)=summary_tables.iter(1,3);
ACGiter(Problem,4)=summary_tables.iter(1,4);
ACGiter(Problem,5)=summary_tables.iter(1,5);

filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');


%%%Problem=22
% Initialize
Problem=22
N = 1000;
r=2;
seed=218;
dimM=20;
dimN=100;
M=10^5;
m=10^3;
global_tol = 1e-5;

[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();

%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;
%%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;

%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam, spa1_hparam, spa2_hparam}; 
name_arr = {'iALM', 'IPL', 'RQP', 'SPA1', 'SPA2'};
framework_arr = {@iALM3, @IAIPAL3,  @penalty3, @sProxALM, @sProxALM};
solver_arr = {@ECG, @ECG, @AIPP, @ECG, @ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);
FVAL(Problem,3)=summary_tables.fval(1,3);
FVAL(Problem,4)=summary_tables.fval(1,4);
FVAL(Problem,5)=summary_tables.fval(1,5);

RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);
RUN(Problem,3)=summary_tables.runtime(1,3);
RUN(Problem,4)=summary_tables.runtime(1,4);
RUN(Problem,5)=summary_tables.runtime(1,5);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);
ACGiter(Problem,3)=summary_tables.iter(1,3);
ACGiter(Problem,4)=summary_tables.iter(1,4);
ACGiter(Problem,5)=summary_tables.iter(1,5);

filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');


%%%Problem=23
% Initialize
Problem=23
N = 1000;
r=5;
seed=218;
dimM=20;
dimN=100;
M=10^5;
m=10^3;
global_tol = 1e-5;

[oracle, hparams] = test_fn_lin_cone_constr_03r(N, r, M, m, seed, dimM, dimN);

% Create the Model object and specify the solver.
ncvx_qsdp = ConstrCompModel(oracle);

% Set the curvatures and the starting point x0.
ncvx_qsdp.x0 = hparams.x0;
ncvx_qsdp.M = hparams.M;
ncvx_qsdp.m = hparams.m;
ncvx_qsdp.K_constr = hparams.K_constr;

% Set the tolerances
ncvx_qsdp.opt_tol = global_tol;
ncvx_qsdp.feas_tol = global_tol;
ncvx_qsdp.time_limit = 3600;

% Add linear constraints
ncvx_qsdp.constr_fn = hparams.constr_fn;
ncvx_qsdp.grad_constr_fn = hparams.grad_constr_fn;
ncvx_qsdp.set_projector = hparams.set_projector;

% Use a relative termination criterion.
ncvx_qsdp.feas_type = 'relative';
ncvx_qsdp.opt_type = 'relative';

% Create some basic hparams.
base_hparam = struct();


%%%% RQP
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';

rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;
%%%%%%%%
iapial_hparam = base_hparam;
iapial_hparam.acg_steptype = 'variable';
iapial_hparam.sigma_min = 0.3;
iapial_hparam.penalty_multiplier = 2;
iapial_hparam.i_reset_multiplier = false;
iapial_hparam.i_reset_prox_center = false;

%%%
% Luo's method.
spa1_hparam = base_hparam;
spa1_hparam.Gamma = 1;
spa2_hparam = base_hparam;
spa2_hparam.Gamma = 10;

% Create the complicated iALM hparams.
ialm_hparam = base_hparam;
ialm_hparam.i_ineq_constr = false;
ialm_hparam.rho0 = hparams.m;
ialm_hparam.L0 = max([hparams.m, hparams.M]);
ialm_hparam.rho_vec = hparams.m_constr_vec;
ialm_hparam.L_vec = hparams.L_constr_vec;
ialm_hparam.B_vec = hparams.K_constr_vec;

% Run a benchmark test and print the summary.
hparam_arr = {ialm_hparam, iapial_hparam, rqp_aipp_hparam, spa1_hparam, spa2_hparam}; 
name_arr = {'iALM', 'IPL', 'RQP', 'SPA1', 'SPA2'};
framework_arr = {@iALM3, @IAIPAL3,  @penalty3, @sProxALM, @sProxALM};
solver_arr = {@ECG, @ECG, @AIPP, @ECG, @ECG};

% Run the test.
[summary_tables, o_mdl] = run_CCM_benchmark(ncvx_qsdp, framework_arr, solver_arr, hparam_arr, name_arr);
o_tbl = [table(dimN, r), summary_tables.all];
disp(o_tbl);


FVAL(Problem,1)=summary_tables.fval(1,1);
FVAL(Problem,2)=summary_tables.fval(1,2);
FVAL(Problem,3)=summary_tables.fval(1,3);
FVAL(Problem,4)=summary_tables.fval(1,4);
FVAL(Problem,5)=summary_tables.fval(1,5);

RUN(Problem,1)=summary_tables.runtime(1,1);
RUN(Problem,2)=summary_tables.runtime(1,2);
RUN(Problem,3)=summary_tables.runtime(1,3);
RUN(Problem,4)=summary_tables.runtime(1,4);
RUN(Problem,5)=summary_tables.runtime(1,5);

ACGiter(Problem,1)=summary_tables.iter(1,1);
ACGiter(Problem,2)=summary_tables.iter(1,2);
ACGiter(Problem,3)=summary_tables.iter(1,3);
ACGiter(Problem,4)=summary_tables.iter(1,4);
ACGiter(Problem,5)=summary_tables.iter(1,5);

filename='BoxConstrained_23Instances_Vec_105.xlsx';
writetable(FVAL,filename,'Sheet',1,'Range','A1');
writetable(RUN,filename,'Sheet',1,'Range','A71');
writetable(ACGiter,filename,'Sheet',1,'Range','A121');




























































