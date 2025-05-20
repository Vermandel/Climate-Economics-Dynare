function [endo_simul] = SSV_sims0(endo_simul,M_,oo_,spfm_exo_simul)
%TRY_EXO_DYN Summary of this function goes here
%   Detailed explanation goes here

persistent id_Z id_gZ id_L id_gL id_rr id_y id_h id_c id_mc id_pi id_w id_SIG id_gSIG id_E id_M id_d id_delthet id_THETA1 id_mu id_tau id_de id_C id_dy id_lb id_pibar;


%if isempty(id_Z)
id_Z   		= strmatch(deblank('Z'),M_.endo_names,'exact');
id_L   		= strmatch(deblank('L'),M_.endo_names,'exact');
id_gZ  		= strmatch(deblank('gZ'),M_.endo_names,'exact');
id_gL  		= strmatch(deblank('gL'),M_.endo_names,'exact');
id_rr  		= strmatch(deblank('r'),M_.endo_names,'exact');
id_y   		= strmatch(deblank('y'),M_.endo_names,'exact');
id_h   		= strmatch(deblank('h'),M_.endo_names,'exact');
id_c   		= strmatch(deblank('c'),M_.endo_names,'exact');
id_mc  		= strmatch(deblank('mc'),M_.endo_names,'exact');
id_pi  		= strmatch(deblank('pi'),M_.endo_names,'exact');
id_w   		= strmatch(deblank('w'),M_.endo_names,'exact');
id_SIG   	= strmatch(deblank('SIG'),M_.endo_names,'exact');
id_gSIG   	= strmatch(deblank('gSIG'),M_.endo_names,'exact');
id_delthet 	= strmatch(deblank('delthet'),M_.endo_names,'exact');
id_THETA1 	= strmatch(deblank('THETA1'),M_.endo_names,'exact');
id_mu    	= strmatch(deblank('mu'),M_.endo_names,'exact');
id_tau   	= strmatch(deblank('tau'),M_.endo_names,'exact');
id_E     	= strmatch(deblank('E'),M_.endo_names,'exact');
id_de     	= strmatch(deblank('de'),M_.endo_names,'exact');
id_M   		= strmatch(deblank('M'),M_.endo_names,'exact');
id_lb   	= strmatch(deblank('lb'),M_.endo_names,'exact');
id_d   		= strmatch(deblank('d'),M_.endo_names,'exact');
id_C   		= strmatch(deblank('C'),M_.endo_names,'exact');
id_dy   	= strmatch(deblank('dy'),M_.endo_names,'exact');
id_pibar   	= strmatch(deblank('pi_bar'),M_.endo_names,'exact');
%end

GZ1   		= M_.params(strmatch(deblank('GZ1'),M_.param_names,'exact'));
LT    		= M_.params(strmatch(deblank('LT'),M_.param_names,'exact'));
lg    		= M_.params(strmatch(deblank('lg'),M_.param_names,'exact'));
beta  		= M_.params(strmatch(deblank('beta'),M_.param_names,'exact'));
sigmaC		= M_.params(strmatch(deblank('sigmaC'),M_.param_names,'exact'));
alpha		= M_.params(strmatch(deblank('alpha'),M_.param_names,'exact'));
sigmaL		= M_.params(strmatch(deblank('sigmaL'),M_.param_names,'exact'));
chi   		= M_.params(strmatch(deblank('chi'),M_.param_names,'exact'));
kappa 		= M_.params(strmatch(deblank('kappa'),M_.param_names,'exact'));
PI    		= M_.params(strmatch(deblank('PI'),M_.param_names,'exact'));
varsigma    = M_.params(strmatch(deblank('varsigma'),M_.param_names,'exact'));
GS1       	= M_.params(strmatch(deblank('GS1'),M_.param_names,'exact'));
xi        	= M_.params(strmatch(deblank('xi'),M_.param_names,'exact'));
gamma     	= M_.params(strmatch(deblank('gamma'),M_.param_names,'exact'));
psi2      	= M_.params(strmatch(deblank('psi2'),M_.param_names,'exact'));
pb        	= M_.params(strmatch(deblank('pb'),M_.param_names,'exact'));
theta2    	= M_.params(strmatch(deblank('theta2'),M_.param_names,'exact'));
deltapb   	= M_.params(strmatch(deblank('deltapb'),M_.param_names,'exact'));
tau0      	= M_.params(strmatch(deblank('tau0'),M_.param_names,'exact'));
delta_M   	= M_.params(strmatch(deblank('delta_M'),M_.param_names,'exact'));
M_1750   	= M_.params(strmatch(deblank('M_1750'),M_.param_names,'exact'));
varphi      = M_.params(strmatch(deblank('varphi'),M_.param_names,'exact'));
nu   	    = M_.params(strmatch(deblank('nu'),M_.param_names,'exact'));
omega  	    = M_.params(strmatch(deblank('omega'),M_.param_names,'exact'));
tf  	    = M_.params(strmatch(deblank('tf'),M_.param_names,'exact'));
rho  	    = M_.params(strmatch(deblank('rho'),M_.param_names,'exact'));
phi_pi  	= M_.params(strmatch(deblank('phi_pi'),M_.param_names,'exact'));
phi_y  	    = M_.params(strmatch(deblank('phi_y'),M_.param_names,'exact'));
lp  	    = M_.params(strmatch(deblank('lp'),M_.param_names,'exact'));
delta_pi_star = M_.params(strmatch(deblank('delta_pi_star'),M_.param_names,'exact'));
pi_regime1  = M_.params(strmatch(deblank('pi_regime1'),M_.param_names,'exact'));
Dc  	    = M_.params(strmatch(deblank('Dc'),M_.param_names,'exact'));

%netPi = @(aa,cc)  (sqrt(aa.^2+2*aa+4*cc+1)+aa-1)*1/2;	% solution of (1+pi)*(pi-a)=c

if nargin > 3 && ~isempty(spfm_exo_simul)
    id_etau  = strmatch(deblank('e_tau'),M_.exo_names,'exact');
end

h0 = endo_simul(id_y,end) ^(1/alpha);

for t=2:(size(endo_simul,2)-1)
    %% TRENDS
    % Z
	endo_simul(id_gZ,t)		= (1-GZ1)*endo_simul(id_gZ,t-1); 
	endo_simul(id_Z,t)      = endo_simul(id_Z,t-1)*(1+endo_simul(id_gZ,t-1)); 
	% L
	endo_simul(id_L,t) 		= endo_simul(id_L,t-1)^(1-lg)*LT^lg;
	endo_simul(id_gL,t)     = endo_simul(id_L,t)/endo_simul(id_L,t-1)-1;
	% SIG
	endo_simul(id_gSIG,t)	= (1-GS1)*endo_simul(id_gSIG,t-1); 
	endo_simul(id_SIG,t)    = endo_simul(id_SIG,t-1)*(1-endo_simul(id_gSIG,t-1)); 
    % THETA
    endo_simul(id_delthet,t) = endo_simul(id_delthet,t-1)*(1-deltapb);
    endo_simul(id_THETA1,t)  = max(pb/1000/theta2*endo_simul(id_delthet,t)*endo_simul(id_SIG,t),0);

	% damages
	d 						= exp(-gamma*(endo_simul(id_M,t-1)-M_1750));
	if ~isempty(id_d)
		endo_simul(id_d,t) = d;
	end
	
    % carbon policy
    tau 					= tau0;
	if ~isempty(id_tau)
		endo_simul(id_tau,t) = tau;
    end
    if ~isempty(id_etau) && t-1 < size(spfm_exo_simul,1) 
		endo_simul(id_tau,t) = endo_simul(id_tau,t)+varphi*spfm_exo_simul(t-1,id_etau);
         tau 				 = tau+varphi*spfm_exo_simul(t-1,id_etau);
    end
	mu	                    = (tau)^(1/(theta2-1));

    if ~isempty(id_mu) 
        endo_simul(id_mu,t) = mu;
    end
	
    if ~isempty(id_pibar)
        endo_simul(id_pibar,t) = PI+( endo_simul(id_pibar,t-1)-PI)*( 1 - delta_pi_star );
        endo_simul(id_pi,t) 	= endo_simul(id_pibar,t);
    else
        endo_simul(id_pi,t) 	= PI;
    end
    
    %
	D 		= d*Dc*endo_simul(id_y,end);

 %   mc 	= (kappa*(1+endo_simul(id_pi,t))*(endo_simul(id_pi,t)-PI)-(1-nu)*beta*(1+endo_simul(id_gZ,t))^(1)*(endo_simul(id_pi,t)-PI)*(1+endo_simul(id_pi,t)) -  (1-varsigma))/varsigma  ;
	c_y                     = 1-endo_simul(id_THETA1,t)*mu^theta2;
    c_temp                  = c_y*d*h0^alpha;
    lb                      = ((c_temp/(1-omega))-omega*(D)/(1-omega))^(-sigmaC);
    mc                      = chi/(lb*(1-omega)^sigmaL)*(c_temp/c_y/d)^((1+sigmaL)/alpha)/(alpha*c_temp/c_y*d)+tf*endo_simul(id_THETA1,t)*mu^theta2+tf*endo_simul(id_THETA1,t)*tau*(1-mu) ;
    if ~isempty(id_mc)
		endo_simul(id_mc,t) = mc;
    end
    
    %endo_simul(id_pi,t) = netPi( PI,  ((1-varsigma)/kappa + varsigma/kappa*mc)/(1 - (1-nu)*beta*(1+endo_simul(id_gZ,t))))  ;


	endo_simul(id_rr,t) 	= (1+endo_simul(id_pi,t))/((1+endo_simul(id_gZ,t))^(-sigmaC)*beta)-1;
	
    % should take into account
%		y 	= d*h^alpha;
%    w                       = d*alpha*y/h*(mc-tf*THETA1*mu^theta2-tf*THETA1*tau*(1-mu));

	h = (d^(1-sigmaC)*(c_y)^-sigmaC*alpha*(1-omega)^sigmaL/chi*(mc-endo_simul(id_THETA1,t-1)*mu^theta2 - endo_simul(id_THETA1,t)*tau*(1-mu)))^(1/(1+sigmaL-alpha*(1-sigmaC)));
	if ~isempty(id_h)
		endo_simul(id_h,t) 		= h;
    end
	
    endo_simul(id_y,t) 		= d*h^alpha;
	endo_simul(id_c,t)		= endo_simul(id_y,t)*c_y;  
    
	endo_simul(id_rr,t) =  (1+endo_simul(id_rr,end))^(1-rho) * (1+endo_simul(id_rr,t-1))^(rho) * ((1+endo_simul(id_pi,t))/((1+PI)))^((1-rho)*phi_pi) * (endo_simul(id_y,t)/endo_simul(id_y,end))^((1-rho)*phi_y) -1;
	
	w      					= alpha*endo_simul(id_y,t)/h*(mc-endo_simul(id_THETA1,t-1)*mu^theta2 - endo_simul(id_THETA1,t)*tau*(1-mu));
	if ~isempty(id_w)
	endo_simul(id_w,t)      = w;
	end
	
    E =  (1-mu)*endo_simul(id_SIG,t)*endo_simul(id_y,t)*endo_simul(id_Z,t)*endo_simul(id_L,t);
    if ~isempty(id_E)
		endo_simul(id_E,t) 	= E;
		endo_simul(id_de,t) =  log((endo_simul(id_E,t)+1e-8)/(endo_simul(id_E,t-1)+1e-8));
    end

	if ~isempty(id_lb)
	   endo_simul(id_lb,t)  		= lb;
	end	
   
	if ~isempty(id_C)
		endo_simul(id_C,t) =  endo_simul(id_c,t)*endo_simul(id_Z,t)*endo_simul(id_L,t);
	end
	
	endo_simul(id_dy,t) =  log((1+endo_simul(id_gZ,t))*(1+endo_simul(id_gL,t))*endo_simul(id_y,t)/endo_simul(id_y,t-1));

    %% CLIMATE BLOCK
	endo_simul(id_M,t) =  M_1750 +(1-delta_M)*(endo_simul(id_M,t-1) - M_1750) + xi*E;

    if ~isreal(endo_simul(:,t))
%           error('lol')
    end
end

end

