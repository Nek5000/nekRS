function TestADIOSRead(verbose)
if (nargin == 0)
    verbose = 0;
end

if ~isnumeric(verbose) || verbose < 0 
    error('TestADIOSRead expects an non-negative integer verbose argument');
end

testHeat(uint8(fix(verbose)))

%
% equalsDbl()
%
function result = equalsDbl(a,b,eps)
result = abs(a-b) < eps;

%
% checkMinMax
%
function checkMinMax(data, expected_min, expected_max, expected_min_index, expected_max_index)

varname=inputname(1);
eps = 0.000001;
[minV, minI] = min(data(:));
[maxV, maxI] = max(data(:));
if (~equalsDbl(expected_min,minV,eps))
    error('ADIOS read test: %s minimum value %s does not match the expected minimum value %s', varname, num2str(minV,8), num2str(expected_min,8));
end
if (~equalsDbl(expected_max,maxV,eps))
    error('ADIOS read test: %s maximum value %s does not match the expected maximum value %s', varname, num2str(maxV,8), num2str(expected_max,8));
end

%
% testHeat
%
function testHeat(verbose)

fprintf('\n************* Test reading heat.bp ************\n');
fprintf('****  Open file...\n');
f=adiosopen('heat.bp','Verbose', verbose);
f.Handlers

fprintf('\n****   Check variables...\n');
assert(size(f.Variables,2)==2, 'Error: Expected number of variables is 2 but got %s', size(f.Variables,2))
assert(strcmp(f.Variables(1).Name, 'T'), 'Error: Second variable is expected to be ''T'' but got %s', f.Variables(2).Name)
assert(strcmp(f.Variables(2).Name, 'dT'), 'Error: First variable is expected to be ''dT'' but got %s', f.Variables(1).Name)

fprintf('\n****   Check attributes...\n');
assert(size(f.Attributes,2)==0, 'Error: Expected number of attributes is 0 but got %s', size(f.Attributes,2))

fprintf('\n****   Read variable data...\n');

T_EntireDataset = adiosread(f,'T');
checkMinMax(T_EntireDataset, 0.0, 200.0, 76, 31);

T_FirstStep = adiosread(f,'T', [1, 1], [12, 8], 1, 1);
checkMinMax(T_FirstStep, 0.0, 200.0, 76, 31);

T_LastStep = adiosread(f,'T', [1, 1], [12, 8], 3, 1);
checkMinMax(T_LastStep, 92.068497597878973693, 102.53127762541851098, 63, 19);



fprintf('\n****   Close file...\n');
adiosclose(f);

