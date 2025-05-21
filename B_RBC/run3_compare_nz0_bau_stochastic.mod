

% Transition versus Business As Usual
% Deterministic (with MIT shocks)
% gauthier@vermandel.fr


@#include 			"model_file.mod"

% set-up the guess for given exo_filtered_ts
[y_guess,oo_,M_]		    = EP_paths_init(oo_,M_,options_,exo_init_ts);

% introduce random shocks
stochastic_dates   = exo_init_ts.dates(2):(exo_init_ts.dates(2)+500);
innovations        = exo_init_ts(stochastic_dates).data;
idS                = find(options_.surprise_shocks);
innovations(:,idS) = randn([length(stochastic_dates) length(idS)])*chol(M_.Sigma_e(idS,idS));
innovations_ts     = exo_init_ts;
innovations_ts    = ds_leftmerge(innovations_ts,dseries(innovations,stochastic_dates(1),M_.exo_names));

% SCENARIO: Business as usual
% Step 1: Compute the deterministic path
[det_bau] 		= EP_deterministic_path(y_guess,exo_init_ts,oo_,M_,options_);
% Step 2: Compute the stochastic path
[sto_bau]   	= EP_stochastic_path(det_bau,innovations_ts,oo_,M_,options_,stochastic_dates);

% SCENARIO: Transition
% Update parameter in implementation of carbon tax
[det_nz0,oo_nz0,M_nz0]    = update_terminal_state('varphi',1,oo_,M_,exo_init_ts,y_guess);
% Step 1: Compute the deterministic path
[det_nz0] 				= EP_deterministic_path(det_nz0,innovations_ts,oo_nz0,M_nz0,options_);
% Step 2: Compute the stochastic path
[sto_nz0]   	= EP_stochastic_path(det_nz0,innovations_ts,oo_nz0,M_nz0,options_,stochastic_dates);




% PLOT FIGURE
nx = 3;ny = 3;varnames=char('y','c','h','d','tau_USD','E','T','M');
nn=size(varnames,1);
figure('Name','Figure 4: Stochastic economy');
for i1 =1:nn
	subplot(nx,ny,i1)
	idx = strmatch(deblank(varnames(i1,:)),M_.endo_names,'exact');
	hold on;
	plotds(varnames(i1,:),sto_nz0,'*-','Color',[64, 156, 14]/255,'LineWidth',1.5,'MarkerIndice',1:50:size(sto_nz0,1))
	plotds(varnames(i1,:),det_nz0,'.--','Color',[156, 40, 14]/255,'LineWidth',1.5,'MarkerIndice',1:50:size(sto_bau,1))
	hold off;
	grid on
	title(['$' M_.endo_names_tex{idx} '$ ' M_.endo_names_long{idx}],'Interpreter','latex')
	xlim([2000 2100])
end
legend('Stochastic Path','Deterministic Path')



% PLOT FIGURE
nx = 3;ny = 3;varnames=char('y','c','h','d','tau_USD','E','T','M');
nn=size(varnames,1);
figure('Name','Figure 5: Stochastic economy');
for i1 =1:nn
	subplot(nx,ny,i1)
	idx = strmatch(deblank(varnames(i1,:)),M_.endo_names,'exact');
	hold on;
	plotds(varnames(i1,:),sto_bau,'Color',[156, 40, 14]/255,'LineWidth',1.5,'MarkerIndice',1:50:size(sto_bau,1))
	plotds(varnames(i1,:),sto_nz0,'*-','Color',[64, 156, 14]/255,'LineWidth',1.5,'MarkerIndice',1:50:size(sto_nz0,1))
	hold off;
	grid on
	title(['$' M_.endo_names_tex{idx} '$ ' M_.endo_names_long{idx}],'Interpreter','latex')
	xlim([2000 2100])
end
legend('Business-As-Usual','Paris-Agreement')



