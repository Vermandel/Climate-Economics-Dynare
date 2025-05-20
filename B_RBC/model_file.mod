

%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------

addpath('../A_functions');

close all;

@#ifndef SIMULATIONS
	@#define SIMULATIONS = 0
@#endif

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

 
var Z   		${Z}$ (long_name='TFP'),  
	L   		${L}$ (long_name='population (billion)'),
	gZ  		${g_{z}}$ (long_name='TFP growth'),  
	SIG 		${\Sigma}$ (long_name='Decoupling rate '),
	THETA1 		${\theta_{1}}$ (long_name='Abatement trend'),
	delthet  	${\delta_{p}}$ (long_name='Abatement depreciation factor'),
    M  			${m}$ (long_name='Carbon stock (GtC) '),
	y 			${\hat{y}}$ (long_name='Detrended GDP'), 
	c 			${\hat{c}}$ (long_name='Detrended C'),  
	r 			${\hat{r}}$ (long_name='Detrended rate'),  
	w 			${\hat{w}}$ (long_name='Detrended wage'), 
	h 			${\hat{h}}$ (long_name='Hours worked'), 
	mu 			${\mu}$ (long_name='Abatement share'), 
	d  			${d_t}$ (long_name='Damage factor'),
	tau 		${\hat{\tau}_t}$ (long_name='Detrended carbon tax')
	E  			${E_t}$ (long_name='Emissions (GtCO2)'),
	T  			${T_t}$ (long_name='Temperatures'),
	tau_USD		${\tau_t}$ (long_name='Carbon tax (\$/C)'),
	s_a 		${\epsilon^A}$ (long_name='Productivity shock');


varexo 	e_a		${\sigma_a}$ 		(long_name='Std productivity')
			; 

parameters   
	M0			$m_{t_{0}}$ 				(long_name='Initial stock of carbon (GtC)'),  
	xi 			$\xi_{m}$ 					(long_name='Marginal atmospheric retention ratio'),  
	LT  		$L_{\infty}$ 				(long_name='Terminal population (billion)'),  
	GZ1  		${delta_{z}}$ 				(long_name='Decay TFP'),  
	Z0 			${Z_{t_{0}}}$ 				(long_name='Initial TFP'),  
	lg			${l_{g}}$ 					(long_name='Population growth'),  
	sigmaC		${\sigma_{c}}$ 				(long_name='Risk aversion'),
	sigmaH		${\sigma_{h}}$ 				(long_name='Labor disutility'),
	beta 		${\beta}$ 					(long_name='Discount factor'),
	chi  		${\chi}$ 					(long_name='Labor disutility'),
	E0  		${e_{t_{0}}}$ 				(long_name='Initial emissions (GtCO$_{2}$)'),  
	Y0  		${y_{t_{0}}}$ 				(long_name='Initial GDP (trillion USD PPP)'),  
	L0  		${l_{t_{0}}}$ 				(long_name='Initial population (billion)'),   
	tau0   		$\tau_{t_{0}}$ 				(long_name='Initial tax'),  
	delthet0  	$\delta_{\theta,t_{0}}$ 	(long_name='Initial decay'),   
	h0  		$h_{t_{0}}$ 				(long_name='Initial hours worked'),    
	mu0   		$\mu_{t_{0}}$ 				(long_name='Initial abatement share'), 
	THETA0   	$\theta_{1,t_{0}}$ 			(long_name='Initial abatement cost-to-gdp'),  
	c0   		$c_{t_{0}}$ 				(long_name='Initial consumption'),  
	y0   		$y_{1,t_{0}}$ 				(long_name='Initial output'),  
	gamma   	$\gamma$ 					(long_name='Climate damage elasticity'),   
	theta2   	$\theta_{2}$ 				(long_name='Abatement cost curvature'),    
	pb			${p_{p}}$ 					(long_name='Cost abatement'),   
	deltapb		${\delta_{pb}}$ 			(long_name='Decay abatement cost'),  
	delta_M		${\delta_{m}}$ 				(long_name='CO$_{2}$ rate of transfer to deep oceans'), 
	gZ0   		$g_{z,t_{0}}$ 				(long_name='Initial TFP growth'),  
	delta_SIG   $g_{\sigma,t_{0}}$ 			(long_name='Initial decoupling rate'), 
	SIG0 
	GS1 
	gZ0_400		${g_{z,t_{0}} \times 400}$ 			(long_name='Initial TFP growth'), 
	gS0_400		${g_{\sigma,t_{0}}}$ 		(long_name='Decay rate decoupling'), 
	GS1_400		${\delta_{\sigma}}$ 		(long_name='Decay rate emission intensity'), 
	rho_a		${\rho_a}$ 					(long_name='Productivity'),
	varphi			${\varphi}$ 				(long_name='Mitigation policy belief'),
	THETA2020
	tau0_USD		${\tau_{t_{0}}}$  (long_name='Initial carbon price (\$/ton)')	
	xi_T
	;

	varexo e_tau;


%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

/* Panel A: Climate Parameters */ 
xi				= 3/11;
gamma			= 2.379e-5;
theta2			= 2.6;
deltapb			= 1-(1-0.017)^(1/4);
delta_M			= 0;

/* Panel B: Economics Parameters */ 
LT				= 10.48;
lg				= 0.025/4;
GZ1				= 0.0072/4;
GS1_400			= 0;

/* Panel C: Initial Conditions */ 
Y0				= 30/4;
E0				= 20.30/4;
THETA2020		= 0.109;
L0				= 4.85;
M0				= 338*2.13-545;
mu0				= 0.0001;
h0				= 1;
rr0_400			= 5;

/* additional parameters */
delthet0		= 1;
y0          	= 1.0;

/* Estimated parameters */ 
rho_a 			= 0.95;
sigmaC			= 1.94787;
sigmaH			= 0.73685;
varphi 			= 0*0.51509;

delta_SIG 			= 0.0033;
gZ0				= 0.0049;
tau0 			= 0.0000;
GS1				= 0;
beta			= 0.995;
xi_T  			= 0.0023;
%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model;
 	%% TREND BLOCK
	[name='productivity trend']
	Z 		= Z(-1)*(1+gZ(-1));
	gZ      = gZ(-1)*(1-GZ1);
	[name='population trend']
	L 		= L(-1)^(1-lg)*LT^lg;	
	[name='emissions trend']
	SIG 	= SIG(-1) * (1-delta_SIG);
	[name='cost of abatement - level']
	THETA1 	= max(pb/1000/theta2*delthet*SIG,0);
	[name='efficiency trend']
	delthet = delthet(-1)*(1-deltapb);
	%% SHOCKS
	[name='shocks']
	s_a = 1-rho_a + rho_a*s_a(-1) + e_a;
	
	%% ENDOGENOUS VARIABLES
	%% HOUSEHOLD
	[name='Interest rate']
	r = 1/beta*((1+gZ(+1))*c(+1)/c)^sigmaC;
	[name='Wage']
	w = chi * h^sigmaH * c^sigmaC;
	%% FIRMS
	[name='technology']
	y 	= s_a*d*h;
	[name='MpH']
	w 	= (1-THETA1*(tau*theta2*(1-mu)+mu^theta2))*d*s_a;
	[name='Optimal abatement']
	mu 	= ((tau0 + varphi*e_tau))^(1/(theta2-1));
	[name='Emissions']	
	E = (1-mu)*SIG*y*Z*L;
	%% GENERAL EQUILIBRIUM
	[name='Ressources Constraint']
	y = c + THETA1*mu^theta2*y;
	%% CLIMATE
	[name='Carbon stock']
	M = (1-delta_M)*M(-1) + xi*E;
	[name='Temperatures']
	T = xi_T*M;
	[name='Optimal abatement']
	tau = tau0 + varphi*e_tau;
	[name='damages']	
	d = exp(-gamma*(M(-1)));
	[name='Tax in USD']	
	tau_USD		= tau*THETA1/(SIG+1e-8)*1000;

end;


%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------


endval;
	% setting initial variables of the simulations
	Z 		= Z0;
	L		= L0;
	gZ		= gZ0;
	M 		= M0;
	delthet	= delthet0;
	THETA1	= THETA0;
	y    	= y0;
	SIG		= SIG0;
end;


steady_state_model;	
	% <---- Baseline Economy -----> %
	%find hours worked = one when no damages
	%d^(1+sigmaC) = chi * h^sigmaH * (h)^sigmaC*;
	chi = 1;
	s_a = 1;
	
	% <----   initial state   ----> %
	SIG0		= E0/((1-mu0)*Y0);
	pb 			= 1000*theta2*THETA2020/(SIG0*(1-delta_SIG)^(4*(2020-1984.75))*(1-deltapb)^(4*(2020-1984.75)));
	THETA0 		= pb/1000/theta2*delthet0*SIG0;
	d0 			= exp(-gamma*M0);
	h0			= ((d0^(1+sigmaC))/chi)^(1/(sigmaH+sigmaC));
	y0 			= d0*h0;
	c0			= y0;
	Z0			= Y0/(L0*y0);
	% <------ terminal state -----> %
	% Trends:
	SIG		= 0;
	THETA1	= 0;
	delthet	= 0;
	tau		= tau0 + varphi*e_tau;
	L		= LT;
	gZ      = 0;
	% Endogenous
	[Z,M]   = get_Z(Z0,gZ0,GZ1,L0,LT,lg,delta_SIG,SIG0,GS1,deltapb,tau0,theta2,pb,chi,sigmaC,xi,M0,gamma,varphi,sigmaH);
	d		= exp(-gamma*(M));
	mu      = ((tau0 + varphi*e_tau))^(1/(theta2-1));
	h		= (d^(1-sigmaC) / chi)^(1/(sigmaH+sigmaC));
	y 		= d*h;
	c 		= y;
	w 		= d;
	r 		= 1/beta;
	T 		= xi_T*M;
	tau_USD	= tau*THETA1/(SIG+1e-8)*1000;
end;

shocks;
	var e_a; stderr 0.007;
end;


%% Some options for simulating the model
options_.initial_guess_path = 'guess_path';							% guess MATLAB file (not compulsory, just speed up estimation)
options_.expectation_window = 100;									% size of expectation window for extended path
options_.forward_path 		= 3000;									% size of simulations from 1984Q4 up to 3000 additional quarters (2777Q3)
options_.ep.Tdrop 			= 0; 									% possibility to drop some initial period before estimating the model
options_.surprise_shocks	= logical(ones(M_.exo_nbr,1));			% define surprises exogenous variables
options_.surprise_shocks(strcmp(M_.exo_names,'e_tau'))=logical(0);	% define announced exogenous variables
options_.dates0 			= dates('1984Q4'); 						% define starting date of simulation y_{0}


% Creating announcement matrix (empty)
exo_init_ts 					= EP_init_shocks(M_,options_);

% If announcements, feel them here
% here expected carbon tax
idtau     						 = strmatch('e_tau',M_.exo_names,'exact');
id_now    						 = find(exo_init_ts.dates==dates('2023Q3')+1);		% initial period of simulation
id_last   						 = find(exo_init_ts.dates==dates('2100Q1'));		% terminal period of simulation
id_taxmax 						 = find(exo_init_ts.dates==dates('2049Q4'));		% when policy is fully implemented
exo_init 						 = exo_init_ts.data;
exo_init(id_now:id_taxmax,idtau) = linspace(0,1,(id_taxmax-id_now)+1);
exo_init(id_taxmax:end,idtau)    = 1;
exo_init(:,idtau) 				 = smoothdata(exo_init(:,idtau),'gaussian',10);
exo_init_ts 					 = dseries(exo_init,exo_init_ts.dates(1),exo_init_ts.name);



%% Save internally the sequence of announced exogenous variables
call_deterministic_exogenous(exo_init_ts,options_);

% Generate a matlab function that update the initial vector y_0
% each time steady state is computed
gen_histval_func(M_);


%% provide the usual checks, conditional on steady state affected by announcements
resid;
steady;
