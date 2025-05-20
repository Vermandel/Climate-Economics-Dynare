function [exo_init_ts] = EP_init_shocks(M_,options_)
%EP_INIT_SHOCKS Summary of this function goes here
%   Detailed explanation goes here
    if ~isfield(options_,'forward_path') || isempty(options_.forward_path) || options_.forward_path==0
        options_.forward_path = 500;
    end
    if ~isfield(options_,'dates0')
        options_.dates0 = dates('1Q1');
    end
    exo_init_ts= dseries(zeros(options_.forward_path,M_.exo_nbr),options_.dates0,M_.exo_names);
end

