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
L0=1;
number=4;

b = 0.001;
nu = 100;
p = 100;
n = 100;
k=20;
s=5;
seed = 301;

m=1/b;
lambda0=1/(2*m);

[f_s,f_n,grad_f_s,prox_psi_n,psi_s,psi_n, grad_psi_s, Lagrangian, params] = test_fn_spca_arnesh2(b, nu, p, n, s, k, seed);

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
pold=zeros(length(params.set_projector(1)),length(params.set_projector(1)));
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
    [znew,v,iter,fistaFailure,Lipfinal]=FISTASPCA4(zold,mu,L,chi,beta,sigma,maxiterFISTA,lambda,c,pold,zold,psi_s,grad_psi_s,prox_psi_n,params,last);
    if fistaFailure==1
        fistaFailedCount=fistaFailedCount+1;
    end
    
    if k-khat==2
        zstorage=zold;
    elseif k-khat==1
        pstorage=pold;
    end
    
    ACG_total_iter=ACG_total_iter+iter;
    
    if psi_s(zold,lambda,c,pold,zold)+psi_n(zold,lambda)<psi_s(znew,lambda,c,pold,zold)+psi_n(znew,lambda)+params.prod_fn(zold-znew,v)
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
        [znew,v,iter,fistaFailure,Lipfinal]=FISTASPCA4(zold,mu,L,chi,beta,sigma,maxiterFISTA,lambda,c,pold,zold,psi_s,grad_psi_s,prox_psi_n,params,last);
        if fistaFailure==1
            fistaFailedCount=fistaFailedCount+1;
        end
        ACG_total_iter=ACG_total_iter+iter;
        ACG_total_iter_failure=ACG_total_iter_failure+iter;
        if  psi_s(zold,lambda,c,pold,zold)+psi_n(zold,lambda)<psi_s(znew,lambda,c,pold,zold)+psi_n(znew,lambda)+params.prod_fn(zold-znew,v)
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
% lambda
% c
end
toc
finalpi=finalz(1:100,1:100);
finalphi=finalz(101:200,1:100);

trace1=trace(finalpi);
trace2=trace(finalphi);
diagonal(:,1)=diag(finalpi);
diagonal(:,2)=diag(finalphi);
funcval=f_s(finalz)+f_n(finalz)
