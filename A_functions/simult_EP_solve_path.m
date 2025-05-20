function [endo_simul,err] = simult_EP_solve_path(exo_simul,endo_simul,Tep_horizon,oo_,M_,options_,dummy_surprises,samplesize)
%SIMULT_EP_SOLVE_PATH Summary of this function goes here
%   Detailed explanation goes here


    % initialize future information (first row is past shocks set to zero by default)
    oo_.exo_simul = zeros(Tep_horizon+1,M_.exo_nbr);
    
    % loop over all innovations
    for t = 2:(samplesize+1)

	    % load states + path
	    oo_.endo_simul = endo_simul(:,(t-1):(t+Tep_horizon-1));
	    % set zeros shock
	    oo_.exo_simul = zeros(size(oo_.endo_simul,2),M_.exo_nbr);
	    % add current structural innovations
	    oo_.exo_simul(2,:) = 0*exo_simul(t,:);
	    % adjust solver
	    options_.periods = size(oo_.endo_simul,2)-(M_.maximum_lead+M_.maximum_lag);
		
        %% solving 
        weight = 1;
		best_weight = 0;
        iter   = 0;
        while iter < 50
            %[t iter err weight]
			iter = iter + 1;
			
			oo_.exo_simul(2,:) = weight*exo_simul(t-1,:) + (1-weight)*zeros(1,M_.exo_nbr);

            if any(~dummy_surprises)
                % add policy announce
                id_det = find(~dummy_surprises);
                oo_.exo_simul(:,id_det) = weight*exo_simul((t-1):(t+Tep_horizon-1),id_det);
            end
			% convex combination of ss
%			oo_.endo_simul(:,[1 Tep_horizon+1]) = weight*endo_simul(:,[(t-1) (t+Tep_horizon-1)]) + (1-weight)*repmat(oo_.steady_state,[1 2]);
			% convec combination with initial state
%			oo_.endo_simul(:,[1]) = weight*endo_simul(:,[(t-1)]) + (1-weight)* endo_simul(:,Tep_horizon+1);
			%oo_.endo_simul(:,[1 Tep_horizon+1]) = weight*endo_simul(:,[(t-1) (t+Tep_horizon-1)]) + (1-weight)*endo_simul0(:,[(t-1) (t+Tep_horizon-1)]);
			
            % call solver
%			oo_ = perfect_foresight_solver_no_global(oo_,M_,options_);
%			oo_ = perfect_foresight_solver_core(M_, options_, oo_);
            [oo_.endo_simul, success, err] = perfect_foresight_solver_core(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, oo_.exo_steady_state, M_, options_);
			
			if success
				endo_simul(:,(t-1):(t+Tep_horizon-1)) = oo_.endo_simul;
				if weight == 1
					break;
				else
					step=max(0.001,(weight-best_weight));
					weight = min([1 best_weight+step+2*(weight-best_weight)]); 
					best_weight = weight;
				end
			else
				oo_.endo_simul = endo_simul(:,(t-1):(t+Tep_horizon-1));
				weight = max([ 0 weight-.5*(weight-best_weight)]); 
			%	error([ int2str(t) ':lol'])
			end	
			
		end
    end
end

