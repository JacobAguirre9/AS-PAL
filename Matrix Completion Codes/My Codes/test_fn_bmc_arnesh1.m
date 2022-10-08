% SPDX-License-Identifier: MIT
% Copyright © 2021 Weiwei "William" Kong

function [f_s,f_n,grad_f_s,prox_psi_n,psi_s,psi_n, grad_psi_s, Lagrangian, params] = test_fn_bmc_arnesh1(data_name, beta, theta, mu, seed)
% Generates the needed functions for the regularized matrix completion problem under box constraints. Requires a data matrix with
% the variable name 'data' that is read from a .mat file with name given by data_name.
%
% Arguments:
%
%   data_name (character vector): Name of the .mat dataset used for the problem.
%
%   beta (double): One of the objective function's hyperparameters.
%
%   theta (double): One of the objective function's hyperparameters.
%
%   mu (double): One of the objective function's hyperparameters.
%
%   seed (int): The number used to seed MATLAB's random number generator.
%
% Returns:
%
%   A pair consisting of an Oracle and a struct. The oracle is first-order oracle underyling the optimization problem and the
%   struct contains the relevant hyperparameters of the problem.
%

% Set the generator.
rng(seed);

% Initialize.
load(data_name, 'data');
[rId, cId] = find(data);
[dim_m, dim_n] = size(data);
k = nnz(data);
P = zeros(dim_m, dim_n);
for i=1:k
    P(rId(i), cId(i)) = 1;
end

% Set the topology (Euclidean).
prod_fn = @(a,b) sum(dot(a, b));
norm_fn = @(a) norm(a, 'fro');

% Other params
rho = beta / theta ^ 2;
data_upper = 5.0;
data_lower = 0.0;

% Constraint
params.constr_fn = @(Z) Z;
params.grad_constr_fn = @(Z, Delta) Delta;
params.K_constr = 1;
params.set_projector = @(Z) max(min(Z, data_upper), data_lower);

% Output params
params.prod_fn = prod_fn;
params.norm_fn = norm_fn;
params.m = 2 * rho * mu;
params.M = 1;
params.x0 = zeros(dim_m, dim_n);


% Parse params
size_data = size(data);
min_rank = min(size_data);

% Computed params
kappa0 = beta / theta;
mu_bar = mu * kappa0;

% Subroutines
kappa = @(x) beta * log(1 + x / theta);
kappa_prime = @(x) beta ./ (theta + x);

s_prox = @(lam,s,L) prox_l1(s, L*lam * mu_bar);

% Create the standard oracle constructs
f_s = @(x,s) 1/2 * norm_fn(P .* (x - data)) ^ 2 + mu * sum(kappa(s) - kappa0 * s);
f_n = @(s) mu_bar * sum(s);
grad_f_s = @(x,s,U,V) P .* (x - data) + U * bsxfun(@times,  mu * (kappa_prime(s) - kappa0), V(:, 1:min_rank)');
prox_psi_n = @(lam,s,L,U,V) U * bsxfun(@times, s_prox(lam,s,L), V(:, 1:min_rank)');

psi_s=@(x,lambda,pen,p,w,s) lambda*(f_s(x,s)-(params.norm_fn(p))^2/(2*pen)+(pen/2)*(params.norm_fn((x+(p/pen))-params.set_projector(x+(p/pen))))^2)+0.5*(norm_fn(x-w))^2;
psi_n=@(lambda,s) lambda*f_n(s);
grad_psi_s=@(x,lambda,pen,p,w,s,U,V) lambda*(grad_f_s(x,s,U,V)+pen*((x+(p/pen))-params.set_projector(x+(p/pen))))+(x-w);
Lagrangian=@(pen,z,p,s) f_s(z,s)+f_n(s)-(params.norm_fn(p))^2/(2*pen)+(pen/2)*(params.norm_fn((z+(p/pen))-params.set_projector(z+(p/pen))))^2;



end
