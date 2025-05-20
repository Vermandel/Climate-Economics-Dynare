function  plotds(y,varargin)
%PLOTDS Summary of this function goes here
%   Detailed explanation goes here

if isdseries(y)
    yds = y;
    vararg_id = 0;
    tosend = 'plot(dseries_to_num(yds),yds.data';
    if deterministic_simuls.size(2)>1
        error('error: too many plot in dseries')
    end
else
    yds = varargin{1};
    vararg_id = 1;
    tosend = ['plot(dseries_to_num(yds),yds.' deblank(y) '.data'];
end


for ix=(1+vararg_id):size(varargin,2)
    tosend = [tosend ', varargin{' int2str(ix) '}'];
end
tosend = [tosend ')'];
eval(tosend)
