function [llk,ee,res,y_,yyy,lnprior,llks] = dsge_llk(theta,dataset_,e_obj,oo_,M_,options_,bayestopt_,varargin)
	%DSGE_INVERT_2O Summary of this function goes here
	%   Detailed explanation goes here
	

	persistent  current_maxllk;
	persistent	best_exo;
	persistent	best_endogenous;
 	persistent	ydet0;   
   
    verbosity = 0;

	oo = oo_;
	oo_.dr.ghx = [];
    %save test; error('dsgellk');
	% initialize and persistent variables
	llk = 0;
    if isempty(current_maxllk)
        
        if verbosity 
              disp('Initialize persistent variables...')
        end

		% initialize
		current_maxllk 	= -inf;
		best_exo 		= zeros(dataset_.size(1),M_.exo_nbr);
        
        if ~isfield(options_,'expectation_window')
	        % horizon of shocks
	        Tep_horizon			= 25;
        else
            Tep_horizon = options_.expectation_window;
        end

        yyy0 = zeros(M_.endo_nbr,Tep_horizon+1,dataset_.size(1));
		ydet0 = zeros(M_.endo_nbr,options_.forward_path);
		
		
        if exist('varargin') && size(varargin,1) > 0
           ee0 = varargin{1};
        end

        if exist('varargin') && size(varargin,1) > 1
           yyy0 = varargin{2};
           disp('yy0 loaded')
        end 

		if exist('ee0') && size(ee0,1) == size(best_exo,1) && size(ee0,2) == M_.exo_nbr
		   best_exo = ee0;
        end
		
    else
        ee0   = best_exo;
        yyy0  = best_endogenous;
    end
	
    if e_obj.update_best_llk == 0 && exist('ee0') && size(ee0,1) == size(best_exo,1) && size(ee0,2) == M_.exo_nbr
        
         if verbosity 
              disp('Load from previous...')
         end
		 
         % replace internal persistent variabme
         if exist('varargin') && size(varargin,1) > 0
            ee0 = varargin{1};
         end
         
		 if exist('varargin') && size(varargin,1) > 2
            yyy0 = varargin{2};
            disp('lol')
         end
		 
    end

	% if theta is not a vector, vectorize it
    if size(theta,2) > 1
        theta = theta';
    end

	% if priors 
    if sum(bayestopt_.pshape)>0 && options_.prior_trunc>0 && e_obj.check_bounds_theta
		bounds = prior_bounds(bayestopt_, options_.prior_trunc)	;	
		theta2 = min(max(bounds.lb,theta),bounds.ub);
		if abs(sum(theta2-theta,'all'))>0
			% then outside prior bounds
			theta=nan(size(theta));
		end
    end

	
	% solve the model/update matrices & steady state
    infos(1)=0;
	if sum(isnan(theta))==0

		% check parameters' interval
		theta2 = min(max(e_obj.theta.lb,theta),e_obj.theta.ub);
		% penalize objective function if bound not restricted
		llk = llk - 10^7*sum(abs(theta2-theta));

		theta=theta2;
		warning off;

		% try update policy functions
		%try

			% copy the structure
			M = M_;
			
			M.params(e_obj.theta.id(e_obj.theta.is_pm)) = theta(e_obj.theta.is_pm);
			M.Sigma_e = zeros(M.exo_nbr);
			%M.Sigma_e(1:length(e_obj.idu.shocks)+1:end)=theta(e_obj.theta.is_sd).^2;
			for ix = e_obj.theta.is_sd
                M.Sigma_e(e_obj.theta.id(ix),e_obj.theta.id(ix)) = theta(ix).^2;
            end
			if isfield(options_,'dynare_vers') && (options_.dynare_vers(2)<6 && options_.dynare_vers(1)==4)
				[oo.dr.ys, M.params,infos(1)] = feval([M.fname '_steadystate2'],oo.dr.ys, zeros(M_.exo_nbr,1), M.params);
			else
				[oo.dr.ys, M.params,infos(1)] = feval([M.fname '.steadystate'],oo.dr.ys, zeros(M_.exo_nbr,1), M.params);
                oo.steady_state = oo.dr.ys;
			end
			
			if ~isfield(options_.ep,'estim')
				% perturbation methods
				[oo.dr,infos(2)] = stochastic_solvers(oo.dr,0,M,options_,oo);
                %[info, oo, options_, M] = stoch_simul(M, options_, oo, {});
			else
				% deterministic simulations
				oo.steady_state=oo.dr.ys;
				% avoid error check for deterministic models
				oo.dr.ghx = 1;
			end

        %catch
       %     disp('dont compute PQ')
		%	infos(end+1)=1;
		%end
	
	else
		infos(end+1)=1;
	end	
	
	% check if solution looks ok
	if (sum(isnan(oo.dr.ys))+sum(isinf(oo.dr.ys))+ sum(abs(imag(oo.dr.ys)))) ~= 0 || sum(oo.dr.ys) > 10^8
		infos(end+1) = 1;
	end

	% if solution has good shape
	if sum(infos)==0
		
		% start likelihood estimation
		N 		= dataset_.size(2);
		T		= dataset_.size(1)-options_.presample;
		T1		= 1+options_.presample;
		% initialize likelihood
		llk		= llk-N*T/2*log(2*pi)-T/2*log(det(M.Sigma_e(e_obj.id.shocks,e_obj.id.shocks)));
		% use the inversion filter to extract the sequence of shocks
		if isfield(options_.ep,'estim')
		%save test;error('lol')
		%tic;[llks,ee,res,y_,stop_dist,yyy,ydet] = IF_EP(dataset_,e_obj,oo,M,options_,ee,yyy,ydet);toc
			% extended path
			[llks,ee,res,y_,stop_dist,yyy,ydet] = IF_EP(dataset_,e_obj,oo,M,options_,ee0,yyy0,ydet0,varargin);		
		elseif options_.order==1 && ~isfield(e_obj,'pw')
			% linearized model
			[llks,ee,res,y_,stop_dist] = IF1(dataset_,e_obj,oo,M,options_,varargin);
		elseif options_.order==1 && isfield(e_obj,'pw')
			% piece-wise linear model
			[llks,ee,res,y_,stop_dist] = IF1_pw(dataset_,e_obj,oo,M,options_,varargin);		
		elseif options_.order==2
			% second order model
			[llks,ee,res,y_,stop_dist] = IF2(dataset_,e_obj,oo,M,options_,varargin);
		elseif options_.order==3
			% third order model
			[llks,ee,res,y_,stop_dist,lambdas] = IF3(dataset_,e_obj,oo,M,options_,ee0);		
		end
        
		% check residuals
		if sum(res) < 0.0001% && stop_dist == 0
			res = 0;
        end
        
		% compute likelihood 
		llk = llk + sum(llks(T1:end)) - (res*options_.penalized_function+stop_dist*options_.penalized_function*10);
		
		% compute priors
        if sum(bayestopt_.pshape)>0
		    lnprior = priordens(theta(e_obj.thet_ids),bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7,bayestopt_.p3,bayestopt_.p4);        
        else
            % likelihood estimation
            lnprior = 0;
        end

		%% check priors
		if lnprior == 0 && sum(bayestopt_.pshape)>0
			disp('Null prior')
		end
		
		%% add priors
		if isnan(lnprior) || isinf(lnprior) || isinf(llk) || isnan(llk)
			llk = -10^12;
		else
			llk = llk + lnprior;
		end
	
	else
		llk = -10^12;
		ee=nan;res=nan;y_=nan;
	end
	
	if current_maxllk <= llk &&  e_obj.update_best_llk
        
        if verbosity
            fprintf('%f < %f : %s \n',current_maxllk,llk,'update it');
        end
        

		
        % save previous
        ve_names        = bayestopt_.name;
        params          = theta;
        current_maxllk  = llk;
        best_filtered   = y_;
		best_exo 		= ee;
		best_endogenous	= yyy;
        eval(['save ' M_.fname '_mle_temp.mat params current_maxllk best_filtered best_endogenous best_exo ve_names e_obj;']);

		if exist('ydet')
			ydet0 = ydet;
		end
		
    else
        if verbosity
            fprintf('%f > %f : %s \n',current_maxllk,llk,'do nothing please');
        end
	end
	
	% minimize the minus
	llk = - llk;
end

