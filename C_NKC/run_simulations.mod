
addpath('../A_functions');

% add additional variables to simulations :
@#define SIMULATIONS = 1
@#include 			"model_file.mod"


% reload estimation
load('run_estimation_filtered_data.mat')
M_ = reload_estimated_parameters(M_,thetamh,theta_names);

steady;
%----------------------------------------------------------------
% Computation
%----------------------------------------------------------------

outsample_dates = dates('2023Q4'):filtered_shock_ts.dates(end);
insample_dates 	= dates('1985Q1'):dates('2023Q3');

%% CHECK SMOOTHED ESTIMATES
%% Filtered sequence
% set-up the guess for given exo_filtered_ts
[deterministic_simuls,oo_,M_]		= EP_paths_init(oo_,M_,options_,exo_init_ts);
% Step 1: Compute the deterministic path
[deterministic_simuls] 				= EP_deterministic_path(deterministic_simuls,exo_init_ts,oo_,M_,options_);
% Step 2: Compute the stochastic path
[stochastic_simuls]   				= EP_stochastic_path(deterministic_simuls,filtered_shock_ts,oo_,M_,options_,insample_dates);
[stochastic_simuls] 				= EP_deterministic_path(stochastic_simuls,filtered_shock_ts,oo_,M_,options_,outsample_dates);

T = dseries_to_num(dataset_);
figure;
for ix=1:dataset_.size(2)
	subplot(1,dataset_.size(2),ix)
	yn = eval(['stochastic_simuls(dataset_.dates).' dataset_.name{ix} '.data']);
	plot(T,dataset_.data(:,ix),T,yn,'--')
	title(dataset_.name{ix})
end
legend('Obs','Smoothed')



% Introduce News shock on realization of carbon tax: Full realization
[endo_nz0,oo_nz0,M_nz0] = update_terminal_state('varphi',1,oo_,M_,filtered_shock_ts,stochastic_simuls);
[det_nz0] 				= EP_deterministic_path(endo_nz0,filtered_shock_ts,oo_nz0,M_nz0,options_,outsample_dates);

% Introduce News shock on realization of carbon tax: No realization
[endo_bau,oo_bau,M_bau] = update_terminal_state('varphi',0,oo_,M_,filtered_shock_ts,stochastic_simuls);
[det_bau] 				= EP_deterministic_path(endo_bau,filtered_shock_ts,oo_bau,M_bau,options_,outsample_dates);



% PLOT FIGURE 3
nx = 3;ny = 3;varnames=char('lny','pi100','rr100','lny_n','output_gap','rrn100','tau_USD','E');
nn=size(varnames,1);
figure('Name','Figure 2: Model-implied projections based on alternative control rates of emissions');
for i1 =1:nn
	subplot(nx,ny,i1)
	idx = strmatch(deblank(varnames(i1,:)),M_.endo_names,'exact');
	hold on;
	plot(dseries_to_num(stochastic_simuls),stochastic_simuls.data(:,idx),'Color',[3, 161, 252]/255),'LineWidth',1.5,'MarkerIndice',1:50:size(stochastic_simuls,1))
	plot(dseries_to_num(endo_nz0),det_nz0.data(:,idx),'^-','Color',[42, 173, 81]/255,'LineWidth',1.5,'MarkerIndice',1:50:size(stochastic_simuls,1))
	plot(dseries_to_num(endo_bau),det_bau.data(:,idx),'-s','Color',[204, 12, 38]/255,'LineWidth',1.5,'MarkerIndice',1:50:size(stochastic_simuls,1))
	hold off;
	grid on
	title(['$' M_.endo_names_tex{idx} '$ ' M_.endo_names_long{idx}],'Interpreter','latex')
	xlim([2022 2100])
end
legend('Estimated tax path','Paris-Agreement','Laissez-faire')




