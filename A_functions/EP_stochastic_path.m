function [updated_path,err] = EP_stochastic_path(endo_path,exo_path,oo_,M_,options_,dates_vec)
    

    %% Read inputs
    if ~isfield(options_,'expectation_window')
	    % horizon of shocks
	    S		= 25;
    else
        S 		= options_.expectation_window;
    end


    %% unpack data
    % First get path over t_0 to t_{0}+T+S
    dates_lag  = (dates_vec(1)-1):(dates_vec(end)+S+1);
    endo_simul = endo_path(dates_lag).data';
    samplesize = dates_vec.length;
    exo_simul  = exo_path(dates_lag).data;
    
    % id announcements
    idA        = find(~options_.surprise_shocks);
    % id surprises
    idS        = find(options_.surprise_shocks);


    % loop over all innovations
    for t = 2:(samplesize+1)
        

        % set new horizon window
		time_vec   = t+(-1:(S+1));

        % load states + path
	    oo_.endo_simul = endo_simul(:,time_vec);

	    % set zeros shock
	    oo_.exo_simul = zeros(size(oo_.endo_simul,2),M_.exo_nbr);

	    % add current structural innovations
	    oo_.exo_simul(2,idS) = exo_simul(t,idS);
        if ~isempty(idA)
            % if any announcement, load them
            oo_.exo_simul(:,idA) = exo_simul(time_vec,idA);
        end
	    % adjust solver
	    options_.periods = size(oo_.endo_simul,2)-(M_.maximum_lead+M_.maximum_lag);
		
        %% solving 
        weight = 1;
		best_weight = 0;
        iter   = 0;
        while iter < 50
            
            %[t iter err weight]
			iter = iter + 1;
			
            % update shocks, adjust difficulty if necessary with weights
			oo_.exo_simul(2,idS) = weight*exo_simul(t,idS) + (1-weight)*zeros(1,length(idS));
            
            % Solve path
            [oo_.endo_simul, success, err] = perfect_foresight_solver_core(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, oo_.exo_steady_state, M_, options_);
			
			if success

				endo_simul(:,time_vec) = oo_.endo_simul;

				if weight == 1

					break;

                else

					step        = max(0.001,(weight-best_weight));
					weight      = min([1 best_weight+step+2*(weight-best_weight)]); 
					best_weight = weight;
                    
                end

            else

				oo_.endo_simul  = endo_simul(:,time_vec);
				weight          = max([ 0 weight-.5*(weight-best_weight)]); 

			end	
			
        end % end while loop

    end % end for loop
    
    % Repack into dseries
    updated_path = endo_path.data ;
    updated_path(find(dates_lag(1)==endo_path.dates):find(dates_lag(end)==endo_path.dates),:) = endo_simul';
    updated_path  = dseries(updated_path,endo_path.dates(1),M_.endo_names);

end