function [left_ds] = ds_leftmerge(left_ds,right_ds)
%DS_LEFTMGE Summary of this function goes here
%   Detailed explanation goes here

    % unpack
    the_data = left_ds.data;

    % left merge
    the_data(find(right_ds.dates(1)==left_ds.dates):find(right_ds.dates(end)==left_ds.dates),:) = right_ds.data;
    
    % repack
    left_ds = dseries(the_data,left_ds.dates(1),left_ds.name);

end

