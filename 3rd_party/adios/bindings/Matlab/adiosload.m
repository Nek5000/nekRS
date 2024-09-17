function adiosload(file, prefix)
%ADIOSLOAD Read all variables in an ADIOS BP file
%
%   ADIOSLOAD is a batch process to load all variables in an ADIOS BP 
%   file. The simplest way is to give an ADIOS BP file name. An optional 
%   'prefix' string can be given to avoid name conflict or name conversion. 
%
%   ADIOSLOAD(FILE)
%      Read all the variables in FILE and use the same variable names.
%
%   ADIOSLOAD(FILE, PREFIX)
%      Read all the variables in FILE and add PREFIX to the variable names.
%
%   See also ADIOSOPEN, ADIOSCLOSE, ADIOS.

%   Copyright 2009 UT-BATTELLE, LLC
%   Date: 2018/10/03
%   Author: Norbert Podhorszki <pnorbert@ornl.gov>

if (~exist('prefix', 'var'))
    prefix = '';
end
fp = adiosopen(file);
for i = 1:length(fp.Variables)
    try
        name{i} = fp.Variables(i).Name;
        data{i} = adiosread(fp, fp.Variables(i).Name);
        assignin('base',[prefix name{i}],data{i});
    catch
        warning(['Skip ... ', fp.Variables(i).Name]);
    end
end
adiosclose(fp);
