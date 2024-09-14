function info = adiosopen(varargin)
%ADIOSOPEN Open an ADIOS BP file and provide information on it.
%   
%   ADIOSOPEN opens an ADIOS BP file and returns a structure that
%   contains the file handler and information on all groups,
%   variables and attributes.
%
%   If you have opened the file with ADIOSOPEN, you can supply the
%   file handler and the group handler to read in a variable/attribute
%   with ADIOSREAD. When finished reading all data, you need to close 
%   the file with ADIOSCLOSE. 
%
%   Note: the adios group is not the same as HDF5 groups. Each variable has
%   a path, which defines a logical hierarchy of the variables within one 
%   adios group. This logical hierarchy is what is similar to HDF5 groups.
%
%   INFO = ADIOSOPEN(FILE) 
%      Open FILE and return an information structure. 
%
%   The returned INFO structure is the following
%
%     Name              File path
%
%     Handlers          Object handlers to pass on to ADIOS functions 
%        FileHandler        uint64 file handler
%        GroupHandler       uint64 IO group object handler
%        ADIOSHandler       uint64 ADIOS object handler
%
%     Variables         Structure array of variables
%           Name            Path of variable
%           Type            Matlab type class of data
%           Dims            Array of dimensions
%           StepsStart      First step's index for this variable in file, always at least 1
%           StepsCount      Number of steps for this variable in file, always at least 1
%           GlobalMin       Global minimum  of the variable (1-by-1 mxArray)
%           GlobalMax       Global maximum of the variable
%           
%     Attribute         Structure array of attributes
%           Name            Path of attribute
%           Type            Matlab type class of data
%           Value           Attribute value
%
%
%   INFO = ADIOSOPEN(FILE, 'Verbose', LEVEL)
%      To get logging from the adiosopen code, set Verbose to 1 or higher.
%      Higher values cause more and more details to be printed.
%
%   See also ADIOSREAD, ADIOSCLOSE, ADIOS.

%   Copyright 2009 UT-BATTELLE, LLC
%   Date: 2018/09/07
%   Author: Norbert Podhorszki <pnorbert@ornl.gov>

%
% Process arguments.
%

checkArgCounts(varargin{:});
[args, msg] = parse_inputs(varargin{:});
if (~isempty(msg))
    error('MATLAB:adiosopen:inputParsing', '%s', msg);
end

if (isnumeric(args.File))
    fn=sprintf('File handler=%lld',args.File);
else
    fn=sprintf('File name=%s',args.File);
end
verbose=sprintf('%d ', args.Verbose);

input = sprintf('adiosopenc:\n  %s \n  Verbose=%s', fn, verbose);
if (args.Verbose > 0) 
    CallArguments = input
end


if (nargout == 1)
    info = adiosopenc(args.File, args.Verbose);
else 
    adiosopenc(args.File, args.Verbose);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:   checkArgCounts %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkArgCounts(varargin)

if (nargin < 1)
    error('MATLAB:adiosopen:notEnoughInputs', ...
          'ADIOSOPEN requires at least one input argument.')
end


if (nargout > 1)
    error('MATLAB:adiosopen:tooManyOutputs', ...
          'ADIOSOPEN requires one or zero output arguments.')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:   parse_inputs   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [args, msg] = parse_inputs(varargin)

args.File    = '';
args.Verbose = 0;  % verbosity, default is off, i.e. 0

msg = '';

% Arg 1: file name
if ischar(varargin{1})
    args.File = varargin{1};
    %
    % Verify existence of filename/dirname
    %
    a = exist(args.File);
    % a is 2 for existing file, 7 for existing directory
    if (a ~= 2 && a ~= 7)
        error('MATLAB:adios:open', ...
              'Couldn''t open file/dir (%d %s).', ...
              a, args.File)
    end
    varargin = {varargin{2:end}};
else
    msg = 'FILE input argument to ADIOSOPEN must be a string ';
    return
end


% Parse optional arguments based on their number.
if (length(varargin) > 0)
    
    paramStrings = {'verbose'};
    
    % For each pair
    for k = 1:2:length(varargin)
        param = lower(varargin{k});
            
        if (~ischar(param))
            msg = 'Parameter name must be a string.';
            return
        end
        
        idx = strmatch(param, paramStrings);
        
        if (isempty(idx))
            msg = sprintf('Unrecognized parameter name "%s".', param);
            return
        elseif (length(idx) > 1)
            msg = sprintf('Ambiguous parameter name "%s".', param);
            return
        end

        switch (paramStrings{idx})
        % VERBOSE
        case 'verbose'
            if (k == length(varargin))
                msg = 'No value specified for Verbose option.';
                return
            end
        
            args.Verbose = varargin{k+1};
            if ((~isnumeric(args.Verbose)) || ...
                (~isempty(find(rem(args.Verbose, 1) ~= 0))))
                
                msg = sprintf('''VERBOSE'' must be an integer.');
                return
            end
            if (args.Verbose < 0)
                msg = sprintf('''VERBOSE'' must be greater or equal to zero.');
                return
            end
        end
    end
end
