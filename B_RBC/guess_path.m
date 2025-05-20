function [endo_simul] = SSV_sims0(endo_simul,M_,oo_,spfm_exo_simul)
%TRY_EXO_DYN Summary of this function goes here
%   Detailed explanation goes here

persistent id_Z id_gZ id_L id_gL id_rr id_y id_h id_c id_mc id_pi id_w id_SIG id_gSIG id_E id_M id_d id_delthet id_THETA1 id_mu id_tau id_de id_C id_dy id_lb id_pibar;


%if isempty(id_Z)
id_Z   		= strmatch(deblank('Z'),M_.endo_names,'exact');
id_L   		= strmatch(deblank('L'),M_.endo_names,'exact');
id_gZ  		= strmatch(deblank('gZ'),M_.endo_names,'exact');
id_rr  		= strmatch(deblank('r'),M_.endo_names,'exact');
id_y   		= strmatch(deblank('y'),M_.endo_names,'exact');
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
id_c   		= strmatch(deblank('c'),M_.endo_names,'exact');
id_h    	= strmatch(deblank('h'),M_.endo_names,'exact');
%end

GZ1   		= M_.params(strmatch(deblank('GZ1'),M_.param_names,'exact'));
LT    		= M_.params(strmatch(deblank('LT'),M_.param_names,'exact'));
lg    		= M_.params(strmatch(deblank('lg'),M_.param_names,'exact'));
beta  		= M_.params(strmatch(deblank('beta'),M_.param_names,'exact'));
sigmaC		= M_.params(strmatch(deblank('sigmaC'),M_.param_names,'exact'));
sigmaH		= M_.params(strmatch(deblank('sigmaH'),M_.param_names,'exact'));
chi   		= M_.params(strmatch(deblank('chi'),M_.param_names,'exact'));
GS1       	= M_.params(strmatch(deblank('GS1'),M_.param_names,'exact'));
xi        	= M_.params(strmatch(deblank('xi'),M_.param_names,'exact'));
gamma     	= M_.params(strmatch(deblank('gamma'),M_.param_names,'exact'));
pb        	= M_.params(strmatch(deblank('pb'),M_.param_names,'exact'));
theta2    	= M_.params(strmatch(deblank('theta2'),M_.param_names,'exact'));
deltapb   	= M_.params(strmatch(deblank('deltapb'),M_.param_names,'exact'));
tau0      	= M_.params(strmatch(deblank('tau0'),M_.param_names,'exact'));
delta_M   	= M_.params(strmatch(deblank('delta_M'),M_.param_names,'exact'));
delta_SIG  	= M_.params(strmatch(deblank('delta_SIG'),M_.param_names,'exact'));
varphi      = M_.params(strmatch(deblank('varphi'),M_.param_names,'exact'));
lp  	    = M_.params(strmatch(deblank('lp'),M_.param_names,'exact'));


if nargin > 3 && ~isempty(spfm_exo_simul)
    id_etau  = strmatch(deblank('e_tau'),M_.exo_names,'exact');
end


for t=2:(size(endo_simul,2)-1)
    %% TRENDS
    % Z
	endo_simul(id_gZ,t)		= (1-GZ1)*endo_simul(id_gZ,t-1); 
	endo_simul(id_Z,t)      = endo_simul(id_Z,t-1)*(1+endo_simul(id_gZ,t-1)); 
	% L
	endo_simul(id_L,t) 		= endo_simul(id_L,t-1)^(1-lg)*LT^lg;
	% SIG
	endo_simul(id_SIG,t)    = endo_simul(id_SIG,t-1)*(1-delta_SIG); 
    % THETA
    endo_simul(id_delthet,t) = endo_simul(id_delthet,t-1)*(1-deltapb);
    endo_simul(id_THETA1,t)  = max(pb/1000/theta2*endo_simul(id_delthet,t)*endo_simul(id_SIG,t),0);

	% damages
	endo_simul(id_d,t)	    = exp(-gamma*(endo_simul(id_M,t-1)));
    endo_simul(id_tau,t)    = tau0+varphi*spfm_exo_simul(t,id_etau);
	endo_simul(id_mu,t)		= endo_simul(id_tau,t)^(1/(theta2-1));
	endo_simul(id_w,t) 	    = (1-endo_simul(id_THETA1,t)*(endo_simul(id_tau,t)*theta2*(1-endo_simul(id_mu,t))+endo_simul(id_mu,t)^theta2))*endo_simul(id_d,t);
	endo_simul(id_y,t)       = (endo_simul(id_w,t)/chi*endo_simul(id_d,t)^sigmaH * (1-endo_simul(id_THETA1,t)*endo_simul(id_mu,t)^theta2)^-sigmaC)^(1/(sigmaC+sigmaH));
	endo_simul(id_h,t)       = endo_simul(id_y,t)/endo_simul(id_d,t);
    endo_simul(id_c,t)       = endo_simul(id_y,t)*(1-endo_simul(id_THETA1,t)*endo_simul(id_mu,t)^theta2);
	endo_simul(id_E,t)      = (1-endo_simul(id_mu,t))*endo_simul(id_SIG,t)*endo_simul(id_y,t)*endo_simul(id_Z,t)*endo_simul(id_L,t);

    %% CLIMATE BLOCK
	endo_simul(id_M,t) =  (1-delta_M)*endo_simul(id_M,t-1) + xi*endo_simul(id_E,t);

    if ~isreal(endo_simul(:,t))
%           error('lol')
    end
end

end

