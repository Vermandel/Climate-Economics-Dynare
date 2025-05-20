function [outputY,err] = EP_solve_path( inputY, inputX, M, options_, oo)
%EP_SOLVE_PATH Summary of this function goes here
%   Detailed explanation goes here

    % use inpputs
    oo.endo_simul  = inputY;
    oo.exo_simul   = inputX;
    
    % set options
    options_.periods = size(oo.endo_simul,2)-(M.maximum_lead+M.maximum_lag);
   
    % solve path
    %oo  = perfect_foresight_solver(M, options_, oo);
    [oo.endo_simul, success, err] = perfect_foresight_solver_core(oo.endo_simul, oo.exo_simul, oo.steady_state, oo.exo_steady_state, M, options_);

    if err < 0.00001
        outputY = oo.endo_simul;
    else
        outputY = inputY;
    end

end

