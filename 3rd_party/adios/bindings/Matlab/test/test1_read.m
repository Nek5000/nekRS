function test1_read(verbose)
if (nargin == 0)
    verbose = 0;
end

if ~isnumeric(verbose) || verbose < 0 
    error('TestADIOSRead expects an non-negative integer verbose argument');
end

testread(uint8(fix(verbose)))

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
if (expected_min_index ~= minI)
    error('ADIOS read test: %s location of minimum value %s does not match the expected location %s', varname, num2str(minI,8), num2str(expected_min_index,8));
end
if (expected_max_index ~= maxI)
    error('ADIOS read test: %s location of maximum value %s does not match the expected location %s', varname, num2str(maxI,8), num2str(expected_max_index,8));
end

%
% testread
%
function testread(verbose)

if (verbose>0)
    fprintf('\n************* Test reading test1.bp ************\n');
    fprintf('****  Open file...\n');
end
f=adiosopen('test1.bp','Verbose', verbose);

if (verbose>0)
    fprintf('\n****   Check variables...\n');
end
assert(size(f.Variables,2)==4, 'Error: Expected number of variables should be 4 but got %s', size(f.Variables,2))
assert(strcmp(f.Variables(1).Name, 'ncols'), 'Error: First variable is expected to be ''ncols'' but got %s', f.Variables(1).Name)
assert(strcmp(f.Variables(2).Name, 'note'), 'Error: First variable is expected to be ''note'' but got %s', f.Variables(2).Name)
assert(strcmp(f.Variables(3).Name, 'nrows'), 'Error: First variable is expected to be ''nrows'' but got %s', f.Variables(3).Name)
assert(strcmp(f.Variables(4).Name, 'temperature2D'), 'Error: Second variable is expected to be ''temperature2D'' but got %s', f.Variables(4).Name)

%fprintf('\n****   Check attributes...\n');
%assert(size(f.Attributes,2)==0, 'Error: Expected number of attributes is 0 but got %s', size(f.Attributes,2))

if (verbose>0)
    fprintf('\n****   Read variable data...\n');
end

nrows = adiosread(f,'nrows');
fprintf('# of rows = %i\n', nrows)

ncols = adiosread(f,'ncols');
fprintf('# of cols = %i\n', ncols)

%note = adiosread(f,'note');
%fprintf('Note: %s\n', note)

temperature2D = adiosread(f,'temperature2D');
checkMinMax(temperature2D, 31.0, 60.0, 1, 30);

fprintf('temperature2d array size = %i\n', size(temperature2D,1)*size(temperature2D,2))
disp(transpose(temperature2D))

fprintf('\n****   Close file...\n');
adiosclose(f);

