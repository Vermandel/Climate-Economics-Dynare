function [deterministic_simuls,err] = EP_deterministic_path(endo_path,exo_path,oo_,M_,options_,dates_vec)
    
    
    if nargin<6 || isempty(dates_vec)
        dates_vec = (endo_path.dates(1)+1):(endo_path.dates(end)-1);
    else
        dates_vec = (dates_vec(1)+1):(dates_vec(end)-1);
    end
  

    %% solve initial path with no shocks
    % get options of simulations
    fulld_vec        = (dates_vec(1)-1):(dates_vec(end)+1);
    Tn               = length(fulld_vec);

    % set shocks (should be zero surprises!)
    oo_.exo_simul    = zeros(Tn,M_.exo_nbr);
    
    
    % add possible announces
    idx              = ~options_.surprise_shocks;
    if any(idx)
        % if any announcement, load them
        Dmin = max(fulld_vec(1),exo_path.dates(1));
        Dmax = min(fulld_vec(end),exo_path.dates(end));
        oo_.exo_simul(find(fulld_vec==Dmin):find(fulld_vec==Dmax),find(idx)) = exo_path(Dmin:Dmax).data(:,find(idx));
    end

    % setup endogenous path
    oo_.endo_simul   = endo_path(fulld_vec).data';

    options_.periods = Tn-(M_.maximum_lead+M_.maximum_lag);
    [endo_simul, success, err] = perfect_foresight_solver_core(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, oo_.exo_steady_state, M_, options_);
    % stack into dseries
    newpath = dseries(endo_simul',fulld_vec(1),M_.endo_names);

    %% merge into input path
    % unpack input
    %thedata=endo_path.data;
    %thedata(find(endo_path.dates==newpath.dates(1)):find(endo_path.dates==newpath.dates(end)),:)=newpath.data;
    % stack into dseries
    %deterministic_simuls = dseries(thedata,endo_path.dates(1),M_.endo_names);
    deterministic_simuls = ds_leftmerge(endo_path,newpath);
end