function [data, attributes] = adiosread(varargin)
%ADIOSREAD Read data from an ADIOS BP file.
%   
%   ADIOSREAD reads data from BP file opened with ADIOSOPEN
%   Provide the structure returned by ADIOSOPEN as the first input argument, and the path to a variable.
%   Assume STRUCT = ADIOSOPEN(filepath).
%   Inspect STRUCT.Variables and STRUCT.Attributes for the list of variables 
%   and attributes available in a file.
%
%   DATA = ADIOSREAD(STRUCT, VARPATH) 
%      Read the entire variable VARPATH from a BP file.
%      STRUCT is the output of ADIOSOPEN.
%      VARPATH is a string to a variable or attribute.
%      If an N-dimensional array variable has multiple steps in the file
%      this function reads all steps and returns an N+1 dimensional array
%      where the last dimension equals the number of steps.
%
%   DATA = ADIOSREAD(STRUCT, INDEX) 
%      Read complete DATA from a BP file with path VARPATH.
%      INDEX points to a variable in the STRUCT.Variables() array. 
%
%   Additional parameters:
%
%   DATA = ADIOSREAD(..., START, COUNT, STEPSTART, STEPCOUNT)
%      Read a portion of a variable. 
%
%      START and COUNT:
%      A slice is defined as two arrays of N integers, where N is the 
%      number of dimensions of the variable, describing the
%      "start" and "count" values. The "start" values start from 1.
%          E.g. [1 5], [10 2] reads the first 10 values in the first dimension
%      and 2 values from the 5th position in the second dimension resulting in
%      a 10-by-2 array. 
%          You can use negative numbers to index from the end of the array
%      as in python. -1 refers to the last element of the array, -2 the one
%      before and so on. 
%          E.g. [-1], [1] reads in the last value of a 1D array. 
%               [1], [-1] reads in the complete 1D array.
%
%      STEPSTART and STEPCOUNT:
%      Similarly, the number of steps from a specific step can be read instead
%      of all data. Steps start from 1. Negative index can be used as well.
%          E.g. -1, 1  will read in the last step from the file
%               n, -1  will read all steps from 'n' to the last one
%
%      Note: for scalars over time, START and COUNT should be [] 
%
%
%   Please read the file adioscopyright.txt for more information.
%
%   See also ADIOSOPEN, ADIOSCLOSE, ADIOS.

%   Copyright 2009 UT-BATTELLE, LLC
%   Date: 2018/09/07
%   Author: Norbert Podhorszki <pnorbert@ornl.gov>

%
% Process arguments.
%

checkArgCounts(varargin{:});
% remember first arg actual name for error prints
Arg1Name = inputname(1);
[args, msg] = parse_inputs(Arg1Name, varargin{:});
if (~isempty(msg))
    error('MATLAB:adiosread:inputParsing', '%s', msg);
end

offsets=sprintf('%d ', args.Starts);
counts=sprintf('%d ', args.Counts);
verbose=sprintf('%d ', args.Verbose);

CallArguments = sprintf('adiosread.m:\n  File name=%s \n  Var=%s\n  Starts=[%s]  Counts=[%s]\n  StepStart=%d  StepCount=%d\n  Verbose=%s', ...
 args.FileName, args.Path, offsets, counts, args.StepStart, args.StepCount, verbose);
if (args.Verbose > 0) 
    CallArguments
end

if (nargout == 1)
    data = adiosreadc(args.File, args.Group, args.Path, args.Starts, args.Counts, args.StepStart, args.StepCount, args.Verbose);
else 
    adiosreadc(args.File, args.Group, args.Path, args.Starts, args.Counts, args.StepStart, args.StepCount, args.Verbose);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:   checkArgCounts %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkArgCounts(varargin)

if (nargin < 2)
    error('MATLAB:adiosread:notEnoughInputs', ...
          'ADIOSREAD requires at least two input arguments.')
end


if (nargout > 1)
    error('MATLAB:adiosread:tooManyOutputs', ...
          'ADIOSREAD requires one or fewer output arguments.')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:   parse_inputs   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [args, msg] = parse_inputs(Arg1Name, varargin)

nargs = nargin - 1; % size of argument list besides the 1st string 

args.File    = uint64(0);   % saved file handler
args.FileName = '';         % file name (for debugging purposes)
args.Group   = uint64(0);   % saved IO group handler
args.ADIOS   = uint64(0);   % saved ADIOS handler
args.Path    = '';          % variable path string
args.VarIndex = -1;         % index of variable in struct.Variables()
args.Starts = int64([]);    % start positions in each dimension for slicing
args.Counts  = int64([]);   % counts in each dimension for slicing
args.StepStart = int64(1);  % starting step
args.StepCount = int64(-1); % number of steps to read
args.Verbose = 0;           % verbosity, default is off, i.e. 0

msg = '';

% Arg 1: struct from ADIOSOPEN
if (isa(varargin{1}, 'struct'))
    try
        infostruct = varargin{1};
        args.File = infostruct.Handlers.FileHandler; % uint64
        args.Group = infostruct.Handlers.GroupHandler; % uint64
        args.ADIOS = infostruct.Handlers.ADIOSHandler; % uint64
    catch
        msg = ['1st argument should be the info struct from ADIOSOPEN'];
        return
    end
else
    msg = ['1st argument should be the info struct from ADIOSOPEN'];
    return
end
args.FileName=infostruct.Name;

% Arg 2: varpath or var index
if (ischar(varargin{2}))
    % VARPATH
    args.Path = varargin{2};
    for k = 1:length(infostruct.Variables)
        if (strcmp(infostruct.Variables(k).Name, args.Path))
            args.VarIndex = k;
            break
        end
    end
    if (args.VarIndex < 0)
        msg = ['2nd argument path does not match any variable names '...
        'in the info struct from ADIOSOPEN'];

    end
elseif (isnumeric(varargin{2}))
    % convert index to int32
    args.VarIndex = int32(varargin{2});
    try 
        % VARPATH
        args.Path = infostruct.Variables(args.VarIndex).Name;
    catch
        msg = ['2nd argument index must be between 1 and number of variables '...
        'in the info struct from ADIOSOPEN'];
        return
    end
else
    msg = ['2nd argument to ADIOSREAD must be a string or a number'];
    return
end

ndim=size(infostruct.Variables(args.VarIndex).Dims,2);

% Arg 3: START array
if (nargs >= 3)
    array = varargin{3};
    if (ndim == 0)
        if (~isempty(array))
            msg = sprintf('3rd argument must be an [] for scalar variables like %s', args.Path);
            return
        end
    else
        if (~isnumeric(array) || isempty(array) || size(array, 1) ~= 1 || ndims(array) ~= 2)
            msg = sprintf('3rd argument must be an 1-by-%u array of integers for variable %s.', ndim, args.Path);
            return
        end
        if(size(array,2) ~= ndim)
           msg = sprintf('3rd argument array size must equal to the dimensions of the variable which is %u in case of variable "%s"', ndim, args.Path);
            return
        end
    end
    args.Starts = int64(fix(array));
end

% Arg 4: COUNT array
if (nargs >= 4)
    array = varargin{4};
    if (ndim == 0)
        if (~isempty(array))
            msg = sprintf('4th argument must be an [] for scalar variables like %s', args.Path);
            return
        end
    else
        if (~isnumeric(array) || isempty(array) || size(array, 1) ~= 1 || ndims(array) ~= 2)
            msg = sprintf('4th argument must be an 1-by-%u array of integers for variable %s.', ndim, args.Path);
            return
        end
        if(size(array,2) ~= ndim)
           msg = sprintf('4th argument array size must equal to the dimensions of the variable which is %u in case of variable "%s"', ndim, args.Path);
            return
        end
    end
    args.Counts = int64(fix(array));
end

% Arg 5: STEPSTART 
if (nargs >= 5)
    value = varargin{5};
    if (~isnumeric(value) || isempty(value) || ndims(value) ~= 2 ||...
       size(value, 1) ~= 1 || size(value, 2) ~= 1 )
        msg = '5th argument must be an 1-by-1 numerical value';
        return
    end
    args.StepStart = int64(fix(value));
    if (args.StepStart == 0 || args.StepStart > infostruct.Variables(args.VarIndex).StepsCount)
       msg = sprintf('5th argument StepStart must be or between 1 and the available Steps (%d for variable "%s)".\nSee %s.Variables(%d).StepsCount', infostruct.Variables(args.VarIndex).StepsCount, args.Path, Arg1Name, args.VarIndex);
       return
    end
    if args.StepStart < 0
        % recalculate negative start here so that we can check stepcount correctly below
        % this calculation is 0..count-1 based here
        while args.StepStart < 0
            args.StepStart = infostruct.Variables(args.VarIndex).StepsCount + args.StepStart;
        end
        % fix back to 1..count base
        args.StepStart = args.StepStart + 1;
    end
end

% Arg 5: STEPCOUNT
if (nargs >= 6)
    value = varargin{6};
    if (~isnumeric(value) || isempty(value) || ndims(value) ~= 2 ||...
       size(value, 1) ~= 1 || size(value, 2) ~= 1 )
        msg = '6th argument must be an 1-by-1 numerical value';
        return
    end
    args.StepCount = int64(fix(value));
    if (args.StepCount == 0)
       msg = '6th argument StepCount cannot be zero';
       return
    end
    if (args.StepStart + args.StepCount - 1  > infostruct.Variables(args.VarIndex).StepsCount)
       msg = sprintf('5th and 6th arguments StepStart and StepCount request steps [%d..%d] beyond the available Steps (%d for variable "%s)".\nSee %s.Variables(%d).StepsCount', args.StepStart, args.StepStart+args.StepCount-1, infostruct.Variables(args.VarIndex).StepsCount, args.Path, Arg1Name, args.VarIndex);
       return
    end
end

