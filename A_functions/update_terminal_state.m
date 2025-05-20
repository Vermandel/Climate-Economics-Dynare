function [endo_path,oo_,M_] = update_terminal_state(pname,pvalue,oo_,M_,exo_path,endo_path)
%TEST Summary of this function goes here
%   Detailed explanation goes here

% update parameter
M_.params(strmatch(pname,M_.param_names,'exact')) = pvalue;

% compute new terminal condition
final_exo = exo_path.data(end,:)';
[oo_.steady_state, M_.params] = eval([ M_.fname  '.steadystate(oo_.steady_state, final_exo, M_.params)']);
yT = oo_.steady_state;
    
% feed terminal value in old path
yT_ts     = dseries(yT',endo_path.dates(end),M_.endo_names);
endo_path = ds_leftmerge(endo_path,yT_ts);

end

