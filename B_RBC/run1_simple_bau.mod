
% Simple Business As Usual Scenario
% Deterministic (no shocks)
% gauthier@vermandel.fr


@#include 			"model_file.mod"

%% Filtered sequence
% set-up the guess for given exo_filtered_ts
[y_guess,oo_,M_]		= EP_paths_init(oo_,M_,options_,exo_init_ts);
% Step 1: Compute the deterministic path
[deterministic_simuls] 				= EP_deterministic_path(y_guess,exo_init_ts,oo_,M_,options_);


% PLOT FIGURE
nx = 1;ny = 4;varnames=char('Z','THETA1','L','SIG');
nn=size(varnames,1);
figure('Name','Figure 1: Exogenous Variables');
for i1 =1:nn
	subplot(nx,ny,i1)
	idx = strmatch(deblank(varnames(i1,:)),M_.endo_names,'exact');
	hold on;
	plotds(varnames(i1,:),deterministic_simuls,'Color',[3, 161, 252]/255,'LineWidth',1.5,'MarkerIndice',1:50:size(deterministic_simuls,1))
	hold off;
	grid on
	title([ M_.endo_names_long{idx} ' - $' M_.endo_names_tex{idx} '$'],'Interpreter','latex')
end


% PLOT FIGURE
nx = 3;ny = 3;varnames=char('y','c','h','d','tau_USD','E','T','M');
nn=size(varnames,1);
figure('Name','Figure 2: Endogenous Variables Business as usual');
for i1 =1:nn
	subplot(nx,ny,i1)
	idx = strmatch(deblank(varnames(i1,:)),M_.endo_names,'exact');
	hold on;
	plotds(varnames(i1,:),deterministic_simuls,'Color',[3, 161, 252]/255,'LineWidth',1.5,'MarkerIndice',1:50:size(deterministic_simuls,1))
	hold off;
	grid on
	title(['$' M_.endo_names_tex{idx} '$ ' M_.endo_names_long{idx}],'Interpreter','latex')
	xlim([2000 2100])
end
legend('Business-as-usual')
