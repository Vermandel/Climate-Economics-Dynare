function [Zend,Mend] = get_Z(Z0,gZ0,GZ1,L0,LT,lg,gS0,SIG0,GS1,deltapb,tau0,theta2,pb,chi,sigmaC,xi,M0,gamma,varphi,sigmaH)
%TRY_EXO_DYN Summary of this function goes here
%   Detailed explanation goes here
    
    e_tau_ts = call_deterministic_exogenous();    
    e_tau = e_tau_ts.e_tau.data;

    T   = length(e_tau);
    gZ  = nan(T,1);
    gZ(1)  = gZ0;

    Z   = nan(T,1);
    Z(1)   = Z0;

    L      = nan(T,1);
    L(1)   = L0;

    SIG    = nan(T,1);
    SIG(1) = SIG0;

    gSIG   = nan(T,1);
    gSIG(1)= gS0;

    delthet = nan(T,1);
    delthet(1)= 1;

    M = nan(T,1);
    M(1) = M0;
    tolx = 1e-6;
    for i = 2:T
        gZ(i)      = (1-GZ1)*gZ(i-1);
        Z(i)       = Z(i-1)*(1+gZ(i-1));
	    L(i)       = L(i-1)^(1-lg)*LT^lg;	
        gSIG(i)    = (1-GS1)*gSIG(i-1);
	    SIG(i) 	   = SIG(i-1)*(1-gSIG(i-1));
	    delthet(i) = delthet(i-1)*(1-deltapb);
	    THETA1     = max(pb/1000/theta2*delthet(i)*SIG(i),0);
        
        tau     = tau0+varphi*e_tau(i);
	    mu		= (tau)^(1/(theta2-1));
	    d		= exp(-gamma*M(i-1));
	    w 	    = (1-THETA1*(tau*theta2*(1-mu)+mu^theta2))*d;
	    y       = (w/chi*d^sigmaH * (1-THETA1*mu^theta2)^-sigmaC)^(1/(sigmaC+sigmaH));
	    h       = y/d;
        c       = y*(1-THETA1*mu^theta2);
        M(i)    = M(i-1) + xi*(1-mu)*SIG(i)*y*Z(i)*L(i);
        if abs(Z(i)-Z(i-1)) < tolx && abs(M(i)-M(i-1)) < tolx
            break;
        end
    end

    Zend = Z(i);
    Mend = M(i);

end

