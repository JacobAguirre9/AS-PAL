% SPDX-License-Identifier: MIT
% Copyright © 2021 Weiwei "William" Kong
% Solve a bounded matrix completion problem instance.
base_hparam = struct();
rqp_aipp_hparam = base_hparam;
rqp_aipp_hparam.aipp_type = 'aipp_v2';
rqp_aipp_hparam.acg_steptype = 'variable';
rqp_aipp_hparam.i_reset_prox_center = false;


Problem=1
% Test generator
beta = 0.5;
theta = 0.5;
mu = 0.5;
seed = 777;
data_name = 'movielens_100k_610u_9724m';
[oracle, hparams] = ...
  test_fn_bmc_01(data_name, beta, theta, mu, seed);

% Create the Model object and specify the solver.
ncvx_bmc = ConstrCompModel(oracle);
ncvx_bmc.solver = @AIPP;
ncvx_bmc.solver_hparams =rqp_aipp_hparam;

ncvx_bmc.feas_type = 'relative';
ncvx_bmc.opt_type = 'relative';

ncvx_bmc.time_limit = 7200;

% Add linear constraints
ncvx_bmc.constr_fn = hparams.constr_fn;
ncvx_bmc.grad_constr_fn = hparams.grad_constr_fn;
ncvx_bmc.set_projector = hparams.set_projector;

% Add penalty framework
ncvx_bmc.K_constr = hparams.K_constr;
ncvx_bmc.opt_tol = 1e-3;
ncvx_bmc.feas_tol = 1e-3;
ncvx_bmc.framework = @penalty500;

% Set the curvatures and the starting point x0.
ncvx_bmc.M = hparams.M;
ncvx_bmc.m = hparams.m;
ncvx_bmc.x0 = hparams.x0;

% Solve the problem.
ncvx_bmc.optimize;
ACGiter(Problem,1)=ncvx_bmc.iter;
RUN(Problem,1)=ncvx_bmc.runtime;

filename='MC_RQP_103_Aug12.xlsx';
writematrix(RUN,filename,'Sheet',1,'Range','A1');
writematrix(ACGiter,filename,'Sheet',1,'Range','A31');


%%%%
Problem=2
% Test generator
beta = 0.5;
theta = 0.5;
mu = 1;
seed = 777;
data_name = 'movielens_100k_610u_9724m';
[oracle, hparams] = ...
  test_fn_bmc_01(data_name, beta, theta, mu, seed);

% Create the Model object and specify the solver.
ncvx_bmc = ConstrCompModel(oracle);
ncvx_bmc.solver = @AIPP;
ncvx_bmc.solver_hparams =rqp_aipp_hparam;

ncvx_bmc.feas_type = 'relative';
ncvx_bmc.opt_type = 'relative';

ncvx_bmc.time_limit = 7200;

% Add linear constraints
ncvx_bmc.constr_fn = hparams.constr_fn;
ncvx_bmc.grad_constr_fn = hparams.grad_constr_fn;
ncvx_bmc.set_projector = hparams.set_projector;

% Add penalty framework
ncvx_bmc.K_constr = hparams.K_constr;
ncvx_bmc.opt_tol = 1e-3;
ncvx_bmc.feas_tol = 1e-3;
ncvx_bmc.framework = @penalty500;

% Set the curvatures and the starting point x0.
ncvx_bmc.M = hparams.M;
ncvx_bmc.m = hparams.m;
ncvx_bmc.x0 = hparams.x0;

% Solve the problem.
ncvx_bmc.optimize;
ACGiter(Problem,1)=ncvx_bmc.iter;
RUN(Problem,1)=ncvx_bmc.runtime;

filename='MC_RQP_103_Aug12.xlsx';
writematrix(RUN,filename,'Sheet',1,'Range','A1');
writematrix(ACGiter,filename,'Sheet',1,'Range','A31');

%%%%
Problem=3
% Test generator
beta = 0.5;
theta = 0.5;
mu = 2;
seed = 777;
data_name = 'movielens_100k_610u_9724m';
[oracle, hparams] = ...
  test_fn_bmc_01(data_name, beta, theta, mu, seed);

% Create the Model object and specify the solver.
ncvx_bmc = ConstrCompModel(oracle);
ncvx_bmc.solver = @AIPP;
ncvx_bmc.solver_hparams =rqp_aipp_hparam;

ncvx_bmc.feas_type = 'relative';
ncvx_bmc.opt_type = 'relative';

ncvx_bmc.time_limit = 7200;

% Add linear constraints
ncvx_bmc.constr_fn = hparams.constr_fn;
ncvx_bmc.grad_constr_fn = hparams.grad_constr_fn;
ncvx_bmc.set_projector = hparams.set_projector;

% Add penalty framework
ncvx_bmc.K_constr = hparams.K_constr;
ncvx_bmc.opt_tol = 1e-3;
ncvx_bmc.feas_tol = 1e-3;
ncvx_bmc.framework = @penalty500;

% Set the curvatures and the starting point x0.
ncvx_bmc.M = hparams.M;
ncvx_bmc.m = hparams.m;
ncvx_bmc.x0 = hparams.x0;

% Solve the problem.
ncvx_bmc.optimize;
ACGiter(Problem,1)=ncvx_bmc.iter;
RUN(Problem,1)=ncvx_bmc.runtime;

filename='MC_RQP_103_Aug12.xlsx';
writematrix(RUN,filename,'Sheet',1,'Range','A1');
writematrix(ACGiter,filename,'Sheet',1,'Range','A31');

%%%%
Problem=4
% Test generator
beta = 0.5;
theta = 1/3;
mu = 0.5;
seed = 777;
data_name = 'movielens_100k_610u_9724m';
[oracle, hparams] = ...
  test_fn_bmc_01(data_name, beta, theta, mu, seed);

% Create the Model object and specify the solver.
ncvx_bmc = ConstrCompModel(oracle);
ncvx_bmc.solver = @AIPP;
ncvx_bmc.solver_hparams =rqp_aipp_hparam;

ncvx_bmc.feas_type = 'relative';
ncvx_bmc.opt_type = 'relative';

ncvx_bmc.time_limit = 7200;

% Add linear constraints
ncvx_bmc.constr_fn = hparams.constr_fn;
ncvx_bmc.grad_constr_fn = hparams.grad_constr_fn;
ncvx_bmc.set_projector = hparams.set_projector;

% Add penalty framework
ncvx_bmc.K_constr = hparams.K_constr;
ncvx_bmc.opt_tol = 1e-3;
ncvx_bmc.feas_tol = 1e-3;
ncvx_bmc.framework = @penalty500;

% Set the curvatures and the starting point x0.
ncvx_bmc.M = hparams.M;
ncvx_bmc.m = hparams.m;
ncvx_bmc.x0 = hparams.x0;

% Solve the problem.
ncvx_bmc.optimize;
ACGiter(Problem,1)=ncvx_bmc.iter;
RUN(Problem,1)=ncvx_bmc.runtime;

filename='MC_RQP_103_Aug12.xlsx';
writematrix(RUN,filename,'Sheet',1,'Range','A1');
writematrix(ACGiter,filename,'Sheet',1,'Range','A31');
%%%%
Problem=5
% Test generator
beta = 0.5;
theta = 1/3;
mu = 1;
seed = 777;
data_name = 'movielens_100k_610u_9724m';
[oracle, hparams] = ...
  test_fn_bmc_01(data_name, beta, theta, mu, seed);

% Create the Model object and specify the solver.
ncvx_bmc = ConstrCompModel(oracle);
ncvx_bmc.solver = @AIPP;
ncvx_bmc.solver_hparams =rqp_aipp_hparam;

ncvx_bmc.feas_type = 'relative';
ncvx_bmc.opt_type = 'relative';

ncvx_bmc.time_limit = 7200;

% Add linear constraints
ncvx_bmc.constr_fn = hparams.constr_fn;
ncvx_bmc.grad_constr_fn = hparams.grad_constr_fn;
ncvx_bmc.set_projector = hparams.set_projector;

% Add penalty framework
ncvx_bmc.K_constr = hparams.K_constr;
ncvx_bmc.opt_tol = 1e-3;
ncvx_bmc.feas_tol = 1e-3;
ncvx_bmc.framework = @penalty500;

% Set the curvatures and the starting point x0.
ncvx_bmc.M = hparams.M;
ncvx_bmc.m = hparams.m;
ncvx_bmc.x0 = hparams.x0;

% Solve the problem.
ncvx_bmc.optimize;
ACGiter(Problem,1)=ncvx_bmc.iter;
RUN(Problem,1)=ncvx_bmc.runtime;

filename='MC_RQP_103_Aug12.xlsx';
writematrix(RUN,filename,'Sheet',1,'Range','A1');
writematrix(ACGiter,filename,'Sheet',1,'Range','A31');

%%%%
Problem=6
% Test generator
beta = 0.5;
theta = 1/3;
mu = 2;
seed = 777;
data_name = 'movielens_100k_610u_9724m';
[oracle, hparams] = ...
  test_fn_bmc_01(data_name, beta, theta, mu, seed);

% Create the Model object and specify the solver.
ncvx_bmc = ConstrCompModel(oracle);
ncvx_bmc.solver = @AIPP;
ncvx_bmc.solver_hparams =rqp_aipp_hparam;

ncvx_bmc.feas_type = 'relative';
ncvx_bmc.opt_type = 'relative';

ncvx_bmc.time_limit = 7200;

% Add linear constraints
ncvx_bmc.constr_fn = hparams.constr_fn;
ncvx_bmc.grad_constr_fn = hparams.grad_constr_fn;
ncvx_bmc.set_projector = hparams.set_projector;

% Add penalty framework
ncvx_bmc.K_constr = hparams.K_constr;
ncvx_bmc.opt_tol = 1e-3;
ncvx_bmc.feas_tol = 1e-3;
ncvx_bmc.framework = @penalty500;

% Set the curvatures and the starting point x0.
ncvx_bmc.M = hparams.M;
ncvx_bmc.m = hparams.m;
ncvx_bmc.x0 = hparams.x0;

% Solve the problem.
ncvx_bmc.optimize;
ACGiter(Problem,1)=ncvx_bmc.iter;
RUN(Problem,1)=ncvx_bmc.runtime;

filename='MC_RQP_103_Aug12.xlsx';
writematrix(RUN,filename,'Sheet',1,'Range','A1');
writematrix(ACGiter,filename,'Sheet',1,'Range','A31');


%%%%
Problem=7
% Test generator
beta = 0.5;
theta = 1/4;
mu = 0.5;
seed = 777;
data_name = 'movielens_100k_610u_9724m';
[oracle, hparams] = ...
  test_fn_bmc_01(data_name, beta, theta, mu, seed);

% Create the Model object and specify the solver.
ncvx_bmc = ConstrCompModel(oracle);
ncvx_bmc.solver = @AIPP;
ncvx_bmc.solver_hparams =rqp_aipp_hparam;

ncvx_bmc.feas_type = 'relative';
ncvx_bmc.opt_type = 'relative';

ncvx_bmc.time_limit = 7200;

% Add linear constraints
ncvx_bmc.constr_fn = hparams.constr_fn;
ncvx_bmc.grad_constr_fn = hparams.grad_constr_fn;
ncvx_bmc.set_projector = hparams.set_projector;

% Add penalty framework
ncvx_bmc.K_constr = hparams.K_constr;
ncvx_bmc.opt_tol = 1e-3;
ncvx_bmc.feas_tol = 1e-3;
ncvx_bmc.framework = @penalty500;

% Set the curvatures and the starting point x0.
ncvx_bmc.M = hparams.M;
ncvx_bmc.m = hparams.m;
ncvx_bmc.x0 = hparams.x0;

% Solve the problem.
ncvx_bmc.optimize;
ACGiter(Problem,1)=ncvx_bmc.iter;
RUN(Problem,1)=ncvx_bmc.runtime;

filename='MC_RQP_103_Aug12.xlsx';
writematrix(RUN,filename,'Sheet',1,'Range','A1');
writematrix(ACGiter,filename,'Sheet',1,'Range','A31');

%%%%
Problem=8
% Test generator
beta = 0.5;
theta = 1/4;
mu = 1;
seed = 777;
data_name = 'movielens_100k_610u_9724m';
[oracle, hparams] = ...
  test_fn_bmc_01(data_name, beta, theta, mu, seed);

% Create the Model object and specify the solver.
ncvx_bmc = ConstrCompModel(oracle);
ncvx_bmc.solver = @AIPP;
ncvx_bmc.solver_hparams =rqp_aipp_hparam;

ncvx_bmc.feas_type = 'relative';
ncvx_bmc.opt_type = 'relative';

ncvx_bmc.time_limit = 7200;

% Add linear constraints
ncvx_bmc.constr_fn = hparams.constr_fn;
ncvx_bmc.grad_constr_fn = hparams.grad_constr_fn;
ncvx_bmc.set_projector = hparams.set_projector;

% Add penalty framework
ncvx_bmc.K_constr = hparams.K_constr;
ncvx_bmc.opt_tol = 1e-3;
ncvx_bmc.feas_tol = 1e-3;
ncvx_bmc.framework = @penalty500;

% Set the curvatures and the starting point x0.
ncvx_bmc.M = hparams.M;
ncvx_bmc.m = hparams.m;
ncvx_bmc.x0 = hparams.x0;

% Solve the problem.
ncvx_bmc.optimize;
ACGiter(Problem,1)=ncvx_bmc.iter;
RUN(Problem,1)=ncvx_bmc.runtime;

filename='MC_RQP_103_Aug12.xlsx';
writematrix(RUN,filename,'Sheet',1,'Range','A1');
writematrix(ACGiter,filename,'Sheet',1,'Range','A31');

%%%%
Problem=9
% Test generator
beta = 0.5;
theta = 0.2;
mu = 0.5;
seed = 777;
data_name = 'movielens_100k_610u_9724m';
[oracle, hparams] = ...
  test_fn_bmc_01(data_name, beta, theta, mu, seed);

% Create the Model object and specify the solver.
ncvx_bmc = ConstrCompModel(oracle);
ncvx_bmc.solver = @AIPP;
ncvx_bmc.solver_hparams =rqp_aipp_hparam;

ncvx_bmc.feas_type = 'relative';
ncvx_bmc.opt_type = 'relative';

ncvx_bmc.time_limit = 7200;

% Add linear constraints
ncvx_bmc.constr_fn = hparams.constr_fn;
ncvx_bmc.grad_constr_fn = hparams.grad_constr_fn;
ncvx_bmc.set_projector = hparams.set_projector;

% Add penalty framework
ncvx_bmc.K_constr = hparams.K_constr;
ncvx_bmc.opt_tol = 1e-3;
ncvx_bmc.feas_tol = 1e-3;
ncvx_bmc.framework = @penalty500;

% Set the curvatures and the starting point x0.
ncvx_bmc.M = hparams.M;
ncvx_bmc.m = hparams.m;
ncvx_bmc.x0 = hparams.x0;

% Solve the problem.
ncvx_bmc.optimize;
ACGiter(Problem,1)=ncvx_bmc.iter;
RUN(Problem,1)=ncvx_bmc.runtime;

filename='MC_RQP_103_Aug12.xlsx';
writematrix(RUN,filename,'Sheet',1,'Range','A1');
writematrix(ACGiter,filename,'Sheet',1,'Range','A31');

%%%%
Problem=10
% Test generator
beta = 0.5;
theta = 0.2;
mu = 1;
seed = 777;
data_name = 'movielens_100k_610u_9724m';
[oracle, hparams] = ...
  test_fn_bmc_01(data_name, beta, theta, mu, seed);

% Create the Model object and specify the solver.
ncvx_bmc = ConstrCompModel(oracle);
ncvx_bmc.solver = @AIPP;
ncvx_bmc.solver_hparams =rqp_aipp_hparam;

ncvx_bmc.feas_type = 'relative';
ncvx_bmc.opt_type = 'relative';

ncvx_bmc.time_limit = 7200;

% Add linear constraints
ncvx_bmc.constr_fn = hparams.constr_fn;
ncvx_bmc.grad_constr_fn = hparams.grad_constr_fn;
ncvx_bmc.set_projector = hparams.set_projector;

% Add penalty framework
ncvx_bmc.K_constr = hparams.K_constr;
ncvx_bmc.opt_tol = 1e-3;
ncvx_bmc.feas_tol = 1e-3;
ncvx_bmc.framework = @penalty500;

% Set the curvatures and the starting point x0.
ncvx_bmc.M = hparams.M;
ncvx_bmc.m = hparams.m;
ncvx_bmc.x0 = hparams.x0;

% Solve the problem.
ncvx_bmc.optimize;
ACGiter(Problem,1)=ncvx_bmc.iter;
RUN(Problem,1)=ncvx_bmc.runtime;

filename='MC_RQP_103_Aug12.xlsx';
writematrix(RUN,filename,'Sheet',1,'Range','A1');
writematrix(ACGiter,filename,'Sheet',1,'Range','A31');

%%%
Problem=11
% Test generator
beta = 0.5;
theta = 1/6;
mu = 0.5;
seed = 777;
data_name = 'movielens_100k_610u_9724m';
[oracle, hparams] = ...
  test_fn_bmc_01(data_name, beta, theta, mu, seed);

% Create the Model object and specify the solver.
ncvx_bmc = ConstrCompModel(oracle);
ncvx_bmc.solver = @AIPP;
ncvx_bmc.solver_hparams =rqp_aipp_hparam;

ncvx_bmc.feas_type = 'relative';
ncvx_bmc.opt_type = 'relative';

ncvx_bmc.time_limit = 7200;

% Add linear constraints
ncvx_bmc.constr_fn = hparams.constr_fn;
ncvx_bmc.grad_constr_fn = hparams.grad_constr_fn;
ncvx_bmc.set_projector = hparams.set_projector;

% Add penalty framework
ncvx_bmc.K_constr = hparams.K_constr;
ncvx_bmc.opt_tol = 1e-3;
ncvx_bmc.feas_tol = 1e-3;
ncvx_bmc.framework = @penalty500;

% Set the curvatures and the starting point x0.
ncvx_bmc.M = hparams.M;
ncvx_bmc.m = hparams.m;
ncvx_bmc.x0 = hparams.x0;

% Solve the problem.
ncvx_bmc.optimize;
ACGiter(Problem,1)=ncvx_bmc.iter;
RUN(Problem,1)=ncvx_bmc.runtime;

filename='MC_RQP_103_Aug12.xlsx';
writematrix(RUN,filename,'Sheet',1,'Range','A1');
writematrix(ACGiter,filename,'Sheet',1,'Range','A31');

Problem=12
% Test generator
beta = 0.5;
theta = 1/6;
mu = 1;
seed = 777;
data_name = 'movielens_100k_610u_9724m';
[oracle, hparams] = ...
  test_fn_bmc_01(data_name, beta, theta, mu, seed);

% Create the Model object and specify the solver.
ncvx_bmc = ConstrCompModel(oracle);
ncvx_bmc.solver = @AIPP;
ncvx_bmc.solver_hparams =rqp_aipp_hparam;

ncvx_bmc.feas_type = 'relative';
ncvx_bmc.opt_type = 'relative';

ncvx_bmc.time_limit = 7200;

% Add linear constraints
ncvx_bmc.constr_fn = hparams.constr_fn;
ncvx_bmc.grad_constr_fn = hparams.grad_constr_fn;
ncvx_bmc.set_projector = hparams.set_projector;

% Add penalty framework
ncvx_bmc.K_constr = hparams.K_constr;
ncvx_bmc.opt_tol = 1e-3;
ncvx_bmc.feas_tol = 1e-3;
ncvx_bmc.framework = @penalty500;

% Set the curvatures and the starting point x0.
ncvx_bmc.M = hparams.M;
ncvx_bmc.m = hparams.m;
ncvx_bmc.x0 = hparams.x0;

% Solve the problem.
ncvx_bmc.optimize;
ACGiter(Problem,1)=ncvx_bmc.iter;
RUN(Problem,1)=ncvx_bmc.runtime;

filename='MC_RQP_103_Aug12.xlsx';
writematrix(RUN,filename,'Sheet',1,'Range','A1');
writematrix(ACGiter,filename,'Sheet',1,'Range','A31');


Problem=13
% Test generator
beta = 0.5;
theta = 1/7;
mu = 0.5;
seed = 777;
data_name = 'movielens_100k_610u_9724m';
[oracle, hparams] = ...
  test_fn_bmc_01(data_name, beta, theta, mu, seed);

% Create the Model object and specify the solver.
ncvx_bmc = ConstrCompModel(oracle);
ncvx_bmc.solver = @AIPP;
ncvx_bmc.solver_hparams =rqp_aipp_hparam;

ncvx_bmc.feas_type = 'relative';
ncvx_bmc.opt_type = 'relative';

ncvx_bmc.time_limit = 7200;

% Add linear constraints
ncvx_bmc.constr_fn = hparams.constr_fn;
ncvx_bmc.grad_constr_fn = hparams.grad_constr_fn;
ncvx_bmc.set_projector = hparams.set_projector;

% Add penalty framework
ncvx_bmc.K_constr = hparams.K_constr;
ncvx_bmc.opt_tol = 1e-3;
ncvx_bmc.feas_tol = 1e-3;
ncvx_bmc.framework = @penalty500;

% Set the curvatures and the starting point x0.
ncvx_bmc.M = hparams.M;
ncvx_bmc.m = hparams.m;
ncvx_bmc.x0 = hparams.x0;

% Solve the problem.
ncvx_bmc.optimize;
ACGiter(Problem,1)=ncvx_bmc.iter;
RUN(Problem,1)=ncvx_bmc.runtime;

filename='MC_RQP_103_Aug12.xlsx';
writematrix(RUN,filename,'Sheet',1,'Range','A1');
writematrix(ACGiter,filename,'Sheet',1,'Range','A31');

Problem=14
% Test generator
beta = 0.5;
theta = 1/7;
mu = 1;
seed = 777;
data_name = 'movielens_100k_610u_9724m';
[oracle, hparams] = ...
  test_fn_bmc_01(data_name, beta, theta, mu, seed);

% Create the Model object and specify the solver.
ncvx_bmc = ConstrCompModel(oracle);
ncvx_bmc.solver = @AIPP;
ncvx_bmc.solver_hparams =rqp_aipp_hparam;

ncvx_bmc.feas_type = 'relative';
ncvx_bmc.opt_type = 'relative';

ncvx_bmc.time_limit = 7200;

% Add linear constraints
ncvx_bmc.constr_fn = hparams.constr_fn;
ncvx_bmc.grad_constr_fn = hparams.grad_constr_fn;
ncvx_bmc.set_projector = hparams.set_projector;

% Add penalty framework
ncvx_bmc.K_constr = hparams.K_constr;
ncvx_bmc.opt_tol = 1e-3;
ncvx_bmc.feas_tol = 1e-3;
ncvx_bmc.framework = @penalty500;

% Set the curvatures and the starting point x0.
ncvx_bmc.M = hparams.M;
ncvx_bmc.m = hparams.m;
ncvx_bmc.x0 = hparams.x0;

% Solve the problem.
ncvx_bmc.optimize;
ACGiter(Problem,1)=ncvx_bmc.iter;
RUN(Problem,1)=ncvx_bmc.runtime;

filename='MC_RQP_103_Aug12.xlsx';
writematrix(RUN,filename,'Sheet',1,'Range','A1');
writematrix(ACGiter,filename,'Sheet',1,'Range','A31');

