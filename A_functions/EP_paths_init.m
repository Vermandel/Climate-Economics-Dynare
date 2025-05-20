function [endo_simul,oo_,M_,yT,y0] = EP_paths_init(oo_,M_,options_,exo_filtered_ts,samplesize)
%EP_PATHS_INIT Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 4 || isempty(exo_filtered_ts)
        % no information about size of simulations, use default ones
        if ~isfield(options_,'forward_path')
            Tepss = 1000;
        else
            Tepss = options_.forward_path;
        end
        if ~isfield(options_,'expectation_window')
	        % horizon of shocks
	        Tep_horizon			= 25;
        else
            Tep_horizon = options_.expectation_window;
        end
        final_exo = zeros(1,M_.exo_nbr);
        init_date = options_.dates0;
        if nargin < 5 || isempty(samplesize)
            samplesize = sum(max(abs(exo_filtered_ts.data(2:end,find(options_.surprise_shocks))),[],2)>0);
        end
        Tsimsize = samplesize+Tep_horizon+Tepss+M_.maximum_lead+M_.maximum_lag;
    else
        % we use 
        final_exo = exo_filtered_ts.data(end,:)';
        init_date = exo_filtered_ts.dates(1);
        Tsimsize  = exo_filtered_ts.size(1);
    end

    if isfield(options_.ep,'Tdrop')
        % take into account possible lags for estimation
        init_date = exo_filtered_ts.dates(1)-options_.ep.Tdrop;
        simsize   = Tsimsize + options_.ep.Tdrop;
    end


    % get terminal state 
    [oo_.steady_state, M_.params] = eval([ M_.fname  '.steadystate(oo_.steady_state, final_exo, M_.params)']);
    y0 = oo_.steady_state;
    yT = oo_.steady_state;
    
    % get initial state if function exists
    if exist(['+' M_.fname '/histval']) > 0
	    % setting initial state
	   y0   	= feval([M_.fname '.histval'], oo_, M_);
    end



    % set initial path
    endo_simul = repmat(y0, [1,Tsimsize]);
    endo_simul(:,end) = yT;

    % load initial guess for path
    if isfield(options_,'initial_guess_path')
        endo_simul = feval(options_.initial_guess_path,endo_simul,M_,oo_,exo_filtered_ts.data);
    end
    
    % set into dseries
    endo_simul = dseries(endo_simul', init_date, M_.endo_names);
end

