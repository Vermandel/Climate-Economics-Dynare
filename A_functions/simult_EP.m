function [endo_simul,endo_simul0,err] = simult_EP(exo_simul,oo_,M_,options_,surprise_shocks,yy,samplesize)
%SIMULT_EP Summary of this function goes here
%   Detailed explanation goes here
    

    options_.noprint = 0;

    if nargin > 4 && ~isempty(surprise_shocks)
        dummy_surprises  = zeros(1,M_.exo_nbr);
        dummy_surprises(surprise_shocks) = 1;
        all_surprises = sort(find(exo_simul(:,surprise_shocks(1))~=0));
        % size of simulation
        try
            samplesize2 = all_surprises(end);
        catch

        end
    else
        dummy_surprises  = ones(1,M_.exo_nbr);
        % size of simulations
        samplesize2 = size(exo_simul,1) ;% size of simulation
    end
    
    if ~isfield(options_,'expectation_window')
	    % horizon of shocks
	    Tep_horizon			= 25;
    else
        Tep_horizon = options_.expectation_window;
    end
    
    if ~isfield(options_,'forward_path')
        Tepss = 1000;
    else
        Tepss = options_.forward_path;
    end

    if ~isfield(options_.ep,'Tdrop')
        options_.ep.Tdrop = 0;
    end   

    if nargin < 7 || isempty(samplesize)
        samplesize = samplesize2;
    end

    oo_.exo_simul = exo_simul(end,:);
    % get terminal state 
    [oo_.steady_state, M_.params] = eval([ M_.fname  '.steadystate(oo_.steady_state, oo_.exo_simul, M_.params)']);
    y0 = oo_.steady_state;

    if exist(['+' M_.fname '/histval']) > 0
	    % setting initial state
	   y0   	= feval([M_.fname '.histval'], oo_, M_);
    end
  
    % set initial path
    endo_simul = repmat(y0, [1,samplesize+Tep_horizon+Tepss+M_.maximum_lead+M_.maximum_lag]);
    endo_simul(:,end) = oo_.steady_state;

    % load initial guess for path
    if isfield(options_,'initial_guess_path')
        endo_simul = feval(options_.initial_guess_path,endo_simul,M_,oo_,exo_simul);
    end
 
    % load initial guess for path
    if nargin>5 && ~isempty(yy)
        endo_simul(:,2:size(yy,2)) =  yy(:,2:size(yy,2));
    end
 
    % solve initial path with no shocks
    Tn               = size(endo_simul,2);
    oo_.exo_simul    = zeros(Tn,M_.exo_nbr);
    if any(~dummy_surprises)
        % if any announcement, load them
        Tn2 = min([size(exo_simul,1) Tn]);
        oo_.exo_simul(1+(1:Tn2),find(~dummy_surprises)) = exo_simul(1:Tn2,find(~dummy_surprises));
    end
    oo_.endo_simul   = endo_simul;
    options_.periods = Tn-(M_.maximum_lead+M_.maximum_lag);
    [oo_.endo_simul, success, err] = perfect_foresight_solver_core(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, oo_.exo_steady_state, M_, options_);

    %perfect_foresight_solver ;
    endo_simul0 = oo_.endo_simul ;
    


    % Get path with no shocks
    endo_simul = oo_.endo_simul(:,options_.ep.Tdrop+(1:(samplesize+1+Tep_horizon)));


    % then solve extended path
    endo_simul = simult_EP_solve_path(exo_simul,endo_simul,Tep_horizon,oo_,M_,options_,dummy_surprises,samplesize);
    
end

