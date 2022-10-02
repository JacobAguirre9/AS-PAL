function [y,u,iter,fistaFailure,L]=FISTASPCA4(x0,mu,L0,chi,beta,sigma,maxiterFISTA,lambda,pen,pold,w,psi_s,grad_psi_s,prox_psi_n,params,last)
iter=0;
x=x0;
y=x0;
A=0;
tau=1;
L=L0;
sumAL=0;
fistaFailure=0;

for j=1:maxiterFISTA
    a=(tau+sqrt(tau^2+4*tau*A*L))/(2*L);
    xtilde=((A*y)+(a*x))/(A+a);
    storegradient=grad_psi_s(xtilde,lambda,pen,pold,w);
    y=prox_psi_n(xtilde-(1/L)*storegradient,lambda, (1/L));
    iter=iter+1;
    
    if psi_s(xtilde,lambda,pen,pold,w)+params.prod_fn(y-xtilde,storegradient)+((1-chi)*L*params.norm_fn(y-xtilde)^2)/(4)+10^-6<psi_s(y,lambda,pen,pold,w)
        success=0;
        L=beta*L;
    else
        success=1;
    end
    while success==0
        a=(tau+sqrt(tau^2+4*tau*A*L))/(2*L);
        xtilde=((A*y)+(a*x))/(A+a);
        storegradient=grad_psi_s(xtilde,lambda,pen,pold,w);
        y=prox_psi_n(xtilde-(1/L)*storegradient,lambda,(1/L));
        iter=iter+1;
        if psi_s(xtilde,lambda,pen,pold,w)+params.prod_fn(y-xtilde,storegradient)+((1-chi)*L*params.norm_fn(y-xtilde)^2)/(4)+10^-6<psi_s(y,lambda,pen,pold,w)
            success=0;
            L=beta*L;
        else
            success=1;
        end
    end
    s=L*(xtilde-y);
    u=grad_psi_s(y,lambda,pen,pold,w)-storegradient+s;
    storegradient=[];
    A=A+a;
    oldtau=tau;
    tau=oldtau+(mu*a);
    x=(1/tau)*((mu*a*y)+(oldtau*x)-(a*s));
    
    
    sumAL=sumAL+A*L*(params.norm_fn(y-xtilde))^2;
    if last==0
        if (params.norm_fn(y-x0))^2<chi*(sumAL)
            fistaFailure=1;
            break
        end
    end
    if last==1
        if (params.norm_fn(y-x0))^2<chi*(A*L*(params.norm_fn(y-xtilde))^2)
            fistaFailure=1;
            break
        end
        
    end
    
    if params.norm_fn(u)<=sigma*params.norm_fn(y-x0)
        break
    else
    end
    
end