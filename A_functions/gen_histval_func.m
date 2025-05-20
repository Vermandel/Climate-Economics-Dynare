function [] = gen_histval_func(M_)
%GEN_HISTVAL_FUNC Summary of this function goes here
%   Detailed explanation goes here

    % Read txt into cell B
    fid = fopen(['+' M_.fname '/driver.m'],'r');
    i = 1;
    tline = fgetl(fid);
    B{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        B{i} = tline;
    end
    fclose(fid);
	

	%% find  '% ENDVAL instructions'
	% in dynare file
	idx = strmatch('% ENDVAL instructions',char(B{:}),'exact');
	% find the next 'steady;' statement
	%idxx = strmatch('steady;',char(B{idx:end}),'exact');
    strn = 'oo_.steady_state';
    idxx = idx +3;
    while strcmp(strn,'oo_.steady_state')
        idxx = idxx +1;
        if size(B{idxx},2) > 15
            strn = B{idxx}(1:16); 
        else
             strn = 'stop';
        end
    end
    
	% extract the relevant content
	func_str = {B{(idx+1):(idxx-1)}};
	% replace oo_.steady_state
	for i1 = 3:size(func_str,2)
		func_str{i1} = strrep(func_str{i1},'oo_.steady_state','ys0_');
	end
	
	% set the name of the function
	func_str{1} = 'function [ys0_] = histval(oo_, M_)';
    func_str{2} = 'ys0_ = oo_.steady_state;';
	% close the function
	func_str{end+1} = 'end';
	
    % Write the Matlab file
    fid = fopen(['+' M_.fname '/histval.m'], 'w');
    for i = 1:numel(func_str)
            fprintf(fid,'%s\n', func_str{i});
    end
	fclose(fid);
end

