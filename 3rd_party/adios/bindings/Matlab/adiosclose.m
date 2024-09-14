function adiosclose(varargin)
%ADIOSCLOSE Close an ADIOS BP file.
%   
%   ADIOSCLOSE closes an ADIOS BP file that was opened by ADIOSOPEN.
%   Provide the structure returned by ADIOSOPEN as input.
%
%   ADIOSCLOSE(STRUCT) 
%      Close file and free internal data structures.
%      STRUCT is the output of ADIOSOPEN.
%
%   ADIOSCLOSE(STRUCT, 'Verbose', LEVEL)
%      To get logging from the adiosclose code, set Verbose to 1 or higher.
%      Higher values cause more and more details to be printed.
%
%   Please read the file adioscopyright.txt for more information.
%
%   See also ADIOSOPEN, ADIOSREAD, ADIOS.

%   Copyright 2009 UT-BATTELLE, LLC
%   Date: 2018/09/07
%   Author: Norbert Podhorszki <pnorbert@ornl.gov>

%
% Process arguments.
%

checkArgCounts(varargin{:});
[args, msg] = parse_inputs(varargin{:});
if (~isempty(msg))
    error('MATLAB:adiosclose:inputParsing', '%s', msg);
end

if (args.Verbose > 0)
    fhs=sprintf('File handler=%lld',args.fh);
    aos=sprintf('ADIOS handler=%lld',args.ao);
    CallArguments = sprintf('adiosclosec:\n  %s\n  %s\n  Verbose=%d\n', fhs, aos, args.Verbose);
end


adiosclosec(args.fh, args.ao, args.Verbose);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:   checkArgCounts %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkArgCounts(varargin)

if (nargin < 1)
    error('MATLAB:adiosopen:notEnoughInputs', ...
          'ADIOSCLOSE requires at least one input argument.')
end


if (nargout > 0)
    error('MATLAB:adiosread:tooManyOutputs', ...
          'ADIOSCLOSE has no output arguments.')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:   parse_inputs   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [args, msg] = parse_inputs(varargin)

arg1  = '';
args.fh=uint64(0);
args.ao=uint64(0);
args.Verbose = 0;  % verbosity, default is off, i.e. 0

msg = '';

% Arg 1: struct returned by ADIOSOPEN
if isstruct(varargin{1})
    arg1 = varargin{1};
    %
    % Verify it is an adios struct from ADIOSOPEN
    %
    fnames = fieldnames(arg1);
    if ( ~strcmp(fnames{1},'Name')        || ...
         ~strcmp(fnames{2},'Handlers') || ...
         ~strcmp(fnames{3},'Variables')   || ...
         ~strcmp(fnames{4},'Attributes') ...
       )
    
        msg = 'STRUCT input argument does not seem to be output of ADIOSOPEN';
        return
    end
    args.fh = uint64(arg1.Handlers.FileHandler);
    args.ao = uint64(arg1.Handlers.ADIOSHandler);
    varargin = {varargin{2:end}};
else
    msg = 'STRUCT input argument must be the output of ADIOSOPEN';
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
