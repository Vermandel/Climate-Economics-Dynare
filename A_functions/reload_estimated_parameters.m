function [M_] = reload_estimated_parameters(M_,theta,theta_names)
%RELOAD_ESTIMATED_PARAMETERS Summary of this function goes here
%   Detailed explanation goes here
    for ix=1:length(theta)
        idx = find(strcmp(deblank(theta_names(ix,:)),M_.param_names));
	    if isempty(idx)
            idx = find(strcmp(deblank(theta_names(ix,:)),M_.exo_names));
            M_.Sigma_e(idx,idx) = theta(ix)^2;
        else
            idx = find(strcmp(deblank(theta_names(ix,:)),M_.param_names));
            M_.params(idx) = theta(ix);
        end
    end
end

