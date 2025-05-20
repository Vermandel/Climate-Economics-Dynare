function [exo_data] = call_deterministic_exogenous(exo_ts,options_)
%CALL_DETERMINISTIC_PATH Summary of this function goes here
%   Detailed explanation goes here
    
    persistent exo_tax;

    if nargin > 0
        % then save tax path
        exo_tax = exo_ts;   
    end
    if nargout > 0
        % then provide saved path
        exo_data = exo_tax;
    end

end

