function [y] = ds_name(y_ds,name)
%DS_NAME Summary of this function goes here
%   Detailed explanation goes here
x='name{1}';
if size(name,1) > 1
    for ix=2:size(name,1) 
       x=[ x ',name{' int2str(ix) '}']; 
    end
end
y = eval(['y_ds{' x '}']);


end

