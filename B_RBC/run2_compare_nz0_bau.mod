
% add additional variables to simulations :
@#include 			"model_file.mod"
% set-up the guess for given exo_filtered_ts
[y_guess,oo_,M_]		    = EP_paths_init(oo_,M_,options_,exo_init_ts);

% SCENARIO: Business as usual
% Step 1: Compute the deterministic path
[y_bau] 				= EP_deterministic_path(y_guess,exo_init_ts,oo_,M_,options_);

% SCENARIO: Transition
% Update parameter in implementation of carbon tax
[y_nz0,oo_nz0,M_nz0]    = update_terminal_state('varphi',1,oo_,M_,exo_init_ts,y_guess);
% Step 1: Compute the deterministic path
[y_nz0] 				= EP_deterministic_path(y_nz0,exo_init_ts,oo_nz0,M_nz0,options_);



% PLOT FIGURE
nx = 3;ny = 3;varnames=char('y','c','h','d','tau_USD','E','T','M');
nn=size(varnames,1);
figure('Name','Figure 2: Model-implied projections based on alternative control rates of emissions');
for i1 =1:nn
	subplot(nx,ny,i1)
	idx = strmatch(deblank(varnames(i1,:)),M_.endo_names,'exact');
	hold on;
	plotds(varnames(i1,:),y_bau,'Color',[156, 40, 14]/255,'LineWidth',1.5,'MarkerIndice',1:50:size(y_bau,1))
	plotds(varnames(i1,:),y_nz0,'*-','Color',[64, 156, 14]/255,'LineWidth',1.5,'MarkerIndice',1:50:size(y_nz0,1))
	hold off;
	grid on
	title(['$' M_.endo_names_tex{idx} '$ ' M_.endo_names_long{idx}],'Interpreter','latex')
	xlim([2000 2100])
end
legend('Business-As-Usual','Paris-Agreement')



