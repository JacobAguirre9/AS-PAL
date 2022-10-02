clc;
clear
close all
sigma=0.1;
rhohat=10^-6;
etahat=10^-6;

mu=0.25;
maxiterFISTA=200000;
maxiter=10000;
chi=0.001;
beta=1.25;
L0=100;
number=75;

N = 1000;
seed = 357;
dimM = 30;
dimN = 100;
M = 10^8;
m = 10^4;
density=0.05;

lambda0=1/(20*m);

%%%% Generate Problem
rng(seed)
NDIMS = 3;
D = sparse(1:dimN, 1:dimN, randi([1, N], 1, dimN), dimN, dimN);
C = sprand(dimM, dimN * dimN, density);
B = sprand(dimN, dimN * dimN, density);
A = sprand(dimM, dimN * dimN, density);
d = rand([dimM, 1]);

% Choose (xi, tau).
[tau, xi, Dfn, Z] = eigen_bisection(M, m, C, D * B);

% Set the topology (Euclidean).
prod_fn = @(a,b) sum(dot(a, b));
norm_fn = @(a) norm(a, 'fro');

% Computing norm of A.
Hp = A * A';
Hp = (Hp + Hp') / 2;
norm_A = sqrt(eigs(Hp, 1, 'la')); % same as lamMax(A'*A)

% Set up helper tensors and operators.
A_tsr = ndSparse(A, [dimM, dimN, dimN]);
At_tsr = permute(A_tsr, NDIMS:-1:1);
B_tsr = ndSparse(B, [dimN, dimN, dimN]);
Bt_tsr = permute(B_tsr, NDIMS:-1:1);
C_tsr = ndSparse(C, [dimM, dimN, dimN]);
Ct_tsr = permute(C_tsr, NDIMS:-1:1);
lin_op = @(Q, x) full(tsr_mult(Q, x, 'primal'));
adj_op = @(Qt, y) sparse(tsr_mult(Qt, y, 'dual'));

% Compute the b vector
E = diag(ones(dimN, 1)) / dimN;
b = lin_op(A_tsr, E);

% Constraint map methods.
params.constr_fn = @(Z) lin_op(A_tsr, Z);
params.grad_constr_fn = @(Z) At_tsr;
params.set_projector = @(Z) b;
params.K_constr = norm_A;

% Basic output params.
A_map = @(Z) lin_op(A_tsr, Z);
params.M = eigs(Dfn(xi, tau) * Z, 1, 'lr');
params.m = -eigs(Dfn(xi, tau) * Z, 1, 'sr');
params.x0 = init_point(dimN, A_map, b, seed);
params.prod_fn = prod_fn;
params.norm_fn = norm_fn;

% Create the Oracle object.
%%William
f_s = @(x) -xi / 2 * norm_fn(D * lin_op(B_tsr, x)) ^ 2 + tau / 2 * norm_fn(lin_op(C_tsr, x) - d) ^ 2;
f_n = @(x) 0;
grad_f_s = @(x) -xi * adj_op(Bt_tsr, (D' * D) * lin_op(B_tsr, x)) + tau * adj_op(Ct_tsr, lin_op(C_tsr, x) - d);
prox_psi_n = @(x, lam) sm_mat_proj(x, 1);
oracle = Oracle(f_s, f_n, grad_f_s, prox_psi_n);

%%%% mine
psi_s=@(x,lambda,pen,p,w) lambda*(f_s(x)+prod_fn(p,A_map(x)-b)+(pen/2)*(norm_fn(A_map(x)-b))^2)+0.5*(norm_fn(x-w))^2;
psi_n=@(x,lambda) lambda*f_n(x);
grad_psi_s=@(x,lambda,pen,p,w)lambda*(grad_f_s(x)+(adj_op(At_tsr,p))+pen*adj_op(At_tsr,(A_map(x)-b)))+(x-w);
Lagrangian=@(pen,z,p) f_s(z)+f_n(z)+params.prod_fn(p,lin_op(A_tsr,z)-b)+(pen/2)*(norm(lin_op(A_tsr,z)-b))^2;

%%%%%%%%%
c1=1;
last=1;


tic
khat=0;
outerIter=0;
cupdate=0;
lower1=1+params.norm_fn(full(grad_f_s(params.x0)));
lower2=1+params.norm_fn(full(params.constr_fn(params.x0)-params.set_projector(1)));

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
Information(20,1)=toc;



%% Generator of an initial point for some of the penalty problems
function x0 = init_point(dimN, A_map, b, seed)

rng(seed);
feasible = true;
while (feasible)
    nVec = 3;
    vMat = rand(dimN, nVec);
    wMat = zeros(dimN, nVec);
    for j=1:nVec
        wMat(:, j) = vMat(:, j) / norm(vMat(:, j));
    end
    lam_unnormed = rand(nVec, 1);
    lam = lam_unnormed / sum(lam_unnormed);
    x0 = zeros(dimN, dimN);
    for j=1:nVec
        x0 = x0 + lam(j) * wMat(:, j) * wMat(:, j)';
    end
    feasible = (norm(A_map(x0) - b, 'fro') <= 1e-6);
end

end




