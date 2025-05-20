function [spfm_exo_simul,options_,Tep_horizon,Tepss] = IF_EP_init_shocks(samplesize,M,options_)
%IF_EP_INIT_SHOCKS Summary of this function goes here
%   Detailed explanation goes here
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

	% initialize future information (first row is past shocks set to zero by default)
	spfm_exo_simul = zeros(2+samplesize+Tep_horizon+Tepss,M.exo_nbr);

end

