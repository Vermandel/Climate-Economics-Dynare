function [T] = dseries_to_num(dataset_)
%DSERIES_TO_NUM Summary of this function goes here
%   Detailed explanation goes here
f = dataset_.freq;
if f==1
    splitter = 'Y';
else
    splitter = 'Q';
end

yyyyqq_t0 = split(char(dataset_.dates(1)),splitter);
yyyyqq_tT = split(char(dataset_.dates(end)),splitter);

T= ((str2num(yyyyqq_t0{1})+(str2num(yyyyqq_t0{2})-1)/f):(1/f):(str2num(yyyyqq_tT{1})+(str2num(yyyyqq_tT{2})-1)/f))';

end

