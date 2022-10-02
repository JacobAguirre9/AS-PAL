clc;
clear
close all
sigma=0.1;
rhohat=10^-5;
etahat=10^-5;

mu=0.25;
maxiterFISTA=200000;
maxiter=10000;
chi=0.001;
beta=1.25;
L0=100;
number=75;

N = 1000;
r=5;
seed = 887;
dimM = 20;
dimN = 100;
M = 10^4;
m = 10^3;
lambda0=20/m;
%%%%Generate Problem
% Initialize.
rng(seed);
D = sparse(1:dimN, 1:dimN, randi([1, N], 1, dimN), dimN, dimN);
C = rand(dimM, dimN);
B = rand(dimN, dimN);
A = rand(dimM, dimN);
d = rand([dimM, 1]);

% Choose (xi, tau).
[tau, xi, ~, ~] = eigen_bisection(M, m, C, D * B);

% Compute the norm of A and other factors.
Hn = B' * (D' * D) * B;
Hp = A' * A;
Hn = (Hn + Hn') / 2;
Hp = (Hp + Hp') / 2;
norm_A = sqrt(eigs(Hp, 1, 'la')); % same as lamMax(A'*A)

% Set the topology (Euclidean)
prod_fn = @(a,b) sum(dot(a, b));
norm_fn = @(a) norm(a, 'fro');

% Compute the b vector;
b = A * (2 * r * rand(dimN, 1) - r);

% Constraint map methods %%%concerning
params.constr_fn = @(z) A * z;
params.grad_constr_fn = @(z) A';
params.set_projector = @(y) b;
params.dual_cone_projector = @(y) y;
params.K_constr = norm_A;

% Basic output params.
params.prod_fn = prod_fn;
params.norm_fn = norm_fn;
params.M = eigs(-xi * Hn + tau * Hp, 1, 'la');
params.m = -eigs(-xi * Hn + tau * Hp, 1, 'sa');
params.x0 = 2 * r * rand(dimN, 1) - r;

% Special params for individual constraints.
params.K_constr_vec = full(sqrt(sum(A .^ 2, 2)));
params.L_constr_vec = zeros(dimM, 1);
params.m_constr_vec = zeros(dimM, 1);

% Oracle construction
f_s = @(x) -xi / 2 * norm_fn(D * B * x) ^ 2 + tau / 2 * norm_fn(C * x - d) ^ 2;
f_n = @(x) 0;
grad_f_s = @(x) -xi * B' * (D' * D) * B * x + tau * C' * (C * x -  d);
prox_psi_n = @(x, lam) box_proj(x, -r, r);
oracle = Oracle(f_s, f_n, grad_f_s, prox_psi_n);

%%%% mine
psi_s=@(x,lambda,pen,p,w) lambda*(f_s(x)+prod_fn(p,A*x-b)+(pen/2)*(norm_fn(A*x-b))^2)+0.5*(norm_fn(x-w))^2;
psi_n=@(x,lambda) lambda*f_n(x);
grad_psi_s=@(x,lambda,pen,p,w)lambda*(grad_f_s(x)+(transpose(A)*p)+pen*transpose(A)*(A*x-b))+(x-w);
Lagrangian=@(pen,z,p) f_s(z)+f_n(z)+params.prod_fn(p,A*z-b)+(pen/2)*(params.norm_fn(A*z-b))^2;

%%%%%%%%%
c1=1;
last=1;

tic
khat=0;
outerIter=0;
cupdate=0;
lower1=1+params.norm_fn(grad_f_s(params.x0));
lower2=1+params.norm_fn(params.constr_fn(params.x0)-params.set_projector(1));

absrhohat=rhohat*lower1;
ACG_total_iter=0;
zold=params.x0;
pold=zeros(length(params.set_projector(1)),1);
C=(2*(1-sigma)^2)/(1-(2*sigma));
c=c1;
lambda=lambda0;
L=L0;
ACG_total_iter_failure=0;
sumlambdaW=0;
sumlambda=0;
ineqFailedCount=0;
fistaFailedCount=0;

for k=1:maxiter
    k
    khat
    outerIter=outerIter+1;
    [znew,v,iter,fistaFailure,Lipfinal]=FISTASDP4(zold,mu,L,chi,beta,sigma,maxiterFISTA,lambda,c,pold,zold,psi_s,grad_psi_s,prox_psi_n,params,last);
    if fistaFailure==1
        fistaFailedCount=fistaFailedCount+1;
    end
    
    if k-khat==2
        zstorage=zold;
    elseif k-khat==1
        pstorage=pold;
    end
    
    ACG_total_iter=ACG_total_iter+iter;
    
    if psi_s(zold,lambda,c,pold,zold)<psi_s(znew,lambda,c,pold,zold)+params.prod_fn(zold-znew,v)
        ineqFailure=1;
        ineqFailedCount=ineqFailedCount+1;
    else
        ineqFailure=0;
    end
    if fistaFailure==1 || ineqFailure==1
        failed=1;
    else
        failed=0;
    end
    
    while fistaFailure==1||ineqFailure==1
        lambda=lambda/2;
        [znew,v,iter,fistaFailure,Lipfinal]=FISTASDP4(zold,mu,L,chi,beta,sigma,maxiterFISTA,lambda,c,pold,zold,psi_s,grad_psi_s,prox_psi_n,params,last);
        if fistaFailure==1
            fistaFailedCount=fistaFailedCount+1;
        end
        ACG_total_iter=ACG_total_iter+iter;
        ACG_total_iter_failure=ACG_total_iter_failure+iter;
        if psi_s(zold,lambda,c,pold,zold)<psi_s(znew,lambda,c,pold,zold)+params.prod_fn(zold-znew,v)
            ineqFailure=1;
            ineqFailedCount=ineqFailedCount+1;
        else
            ineqFailure=0;
        end
    end
    L=Lipfinal
    ACG_failure_percentage=(ACG_total_iter_failure/ACG_total_iter)*100;
    pnew=pold+c*(params.constr_fn(znew)-params.set_projector(1));
    
    w=(v+zold-znew)/(lambda);
    normw=(params.norm_fn(w))/lower1
    normfeas=(params.norm_fn(params.constr_fn(znew)-params.set_projector(1)))/(lower2)
    if normw<=rhohat && normfeas<=etahat
        finalz=znew;
        break
    end
    
    if k>=khat+2
        sumlambdaW=sumlambdaW+(lambda*(params.norm_fn(w))^2);
        sumlambda=sumlambda+lambda;
    end
    
    if k>=khat+2 && (Lagrangian(c,zstorage,pstorage)-Lagrangian(c,znew,pnew)-(params.norm_fn(pnew))^2/(2*c))<=max((sumlambdaW)/(2*C),sumlambda*((absrhohat^2)/(2*C)))
        numberofoutercycle=k;
        c=2*c;
        cupdate=cupdate+1;
        khat=k;
        sumlambdaW=0;
        sumlambda=0;
    else
        c=c;
    end
    zold=znew;
    pold=pnew;
    if iter<=number && failed==0
        lambda=2*lambda;
    end
end
toc

Information(1,1)=0;
Information(2,1)=0;
Information(3,1)=0;
Information(4,1)=0;
Information(5,1)=0;
Information(6,1)=0;
Information(7,1)=lambda;
Information(8,1)=ACG_total_iter;
Information(9,1)=ACG_total_iter_failure;
Information(10,1)=ACG_failure_percentage;
Information(11,1)=outerIter;
Information(12,1)=L;
Information(13,1)=cupdate;
Information(14,1)=c;
Information(15,1)=fistaFailedCount;
Information(16,1)=ineqFailedCount;


Information(1,1)=seed;
Information(2,1)=dimM;
Information(3,1)=dimN;
Information(4,1)=M;
Information(5,1)=m;
Information(6,1)=lambda0;
Information(17,1)=Information(11,1)/(Information(13,1)+1);
Information(18,1)=Information(8,1)/(Information(11,1));
Information(19,1)=f_s(finalz);
Information(20,1)=r;
Information(21,1)=toc;
