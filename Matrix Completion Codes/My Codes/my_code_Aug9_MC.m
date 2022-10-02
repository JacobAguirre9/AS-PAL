clc;
clear
close all
format long
sigma=0.1;
rhohat=10^-3;
etahat=10^-3;

mu=0.25;
maxiterFISTA=200000;
maxiter=10000;
chi=0.001;
beta=1.25;
L0=1;
number=4;

betastar=0.5;
thetastar=1/7;
mustar=1;

m=(2*betastar*mustar)/(thetastar^2);
%lambda0=1/(2*m);
lambda0=10/m;
%lambda0=1;

[f_s,f_n,grad_f_s,prox_psi_n,psi_s,psi_n, grad_psi_s, Lagrangian, params] = test_fn_bmc_arnesh1('movielens_100k_610u_9724m',0.5,1/7,1,777);

%%%%%%%%%
c1=500;
last=1;

tic
khat=0;
outerIter=0;
cupdate=0;
[Uold,Sold,Vold]=svd(params.x0,'econ');
sold=diag(Sold);
lower1=1+params.norm_fn(grad_f_s(params.x0,sold,Uold,Vold));
lower2=1+params.norm_fn(params.x0-params.set_projector(params.x0));
Uold=[];
Vold=[];
Sold=[];

absrhohat=rhohat*lower1;
ACG_total_iter=0;
zold=params.x0;
pold=zeros(610,9724);
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
    [znew,v,iter,fistaFailure,Lipfinal]=FISTAMC1(zold,mu,L,chi,beta,sigma,maxiterFISTA,lambda,c,pold,zold,psi_s,grad_psi_s,prox_psi_n,params,last);
    iter
    if fistaFailure==1
        fistaFailedCount=fistaFailedCount+1;
    end
    
    if k-khat==2
        zstorage=zold;
        s_storage=sold;
    elseif k-khat==1
        pstorage=pold;
    end
    
    ACG_total_iter=ACG_total_iter+iter;
    
    snew=diag(diag(svd(znew,'econ')));
    if psi_s(zold,lambda,c,pold,zold,sold)+psi_n(lambda,sold)<psi_s(znew,lambda,c,pold,zold,snew)+psi_n(lambda,snew)+params.prod_fn(zold-znew,v)
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
        [znew,v,iter,fistaFailure,Lipfinal]=FISTAMC1(zold,mu,L,chi,beta,sigma,maxiterFISTA,lambda,c,pold,zold,psi_s,grad_psi_s,prox_psi_n,params,last);
        if fistaFailure==1
            fistaFailedCount=fistaFailedCount+1;
        end
        ACG_total_iter=ACG_total_iter+iter;
        ACG_total_iter_failure=ACG_total_iter_failure+iter;
       
        snew=diag(diag(svd(znew,'econ')));
        if  psi_s(zold,lambda,c,pold,zold,sold)+psi_n(lambda,sold)<psi_s(znew,lambda,c,pold,zold,snew)+psi_n(lambda,snew)+params.prod_fn(zold-znew,v)
            ineqFailure=1;
            ineqFailedCount=ineqFailedCount+1;
        else
            ineqFailure=0;
        end
    end
    L=Lipfinal
    ACG_failure_percentage=(ACG_total_iter_failure/ACG_total_iter)*100;
    pnew=pold+c*(znew-params.set_projector(znew+(pold/c)));
    
    w=(v+zold-znew)/(lambda);
    normw=(params.norm_fn(w))/lower1
    normfeas=(params.norm_fn((znew+pold/c)-params.set_projector(znew+pold/c)))/(lower2)
    if normw<=rhohat && normfeas<=etahat
        finalz=znew;
        finals=diag(diag(svd(znew,'econ')));
        funcvalue=f_s(finalz,finals)+f_n(finals);
        break
    end
    
    if k>=khat+2
        sumlambdaW=sumlambdaW+(lambda*(params.norm_fn(w))^2);
        sumlambda=sumlambda+lambda;
    end
    if k>=khat+2 && (Lagrangian(c,zstorage,pstorage,s_storage)-Lagrangian(c,znew,pnew,snew)-(params.norm_fn(pnew))^2/(2*c))<=max((sumlambdaW)/(2*C),sumlambda*((absrhohat^2)/(2*C)))
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
    sold=snew;
    pold=pnew;
   if iter<=number && failed==0
        lambda=2*lambda;
   end
lambda
% c
end
toc