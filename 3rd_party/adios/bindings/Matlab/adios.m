function adios()
%ADIOS Reader for the ADIOS BP file format
%   
%   ADIOS is a componentized I/O layer for high performance output combined
%   with an easy to use interface and a set of different data writing methods
%   so that everyone can choose the best method for a particular system.
%
%   ADIOSOPEN     returns a structure with information on an ADIOS BP File
%                 variables and attributes.
%   ADIOSREAD     reads in a variable from the file.
%                 It expects the info structure returned by ADIOSOPEN.
%   ADIOSREADATTR reads in an attribute from the file.
%                 It expects the info structure returned by ADIOSOPEN.
%   ADIOSCLOSE    closes the file.
%
%   Organization of an ADIOS BP file
%   --------------------------------
%   An ADIOS BP file contains a set of variables and attributes.
%   Each variable in the group has a path, which defines a logical 
%   hierarchy of the variables within the file. 
%
%   Time dimension of a variable
%   ----------------------------
%   Variables can be written several times from a program, if they have a time
%   dimension. The reader exposes the variables with an extra dimension, i.e.
%   a 2D variable written over time is seen as a 3D variable. In MATLAB, the
%   extra dimension is the last dimension (the slowest changing dimension).
%   Since the reader allows reading an arbitrary slice of a variable, data for
%   one timestep can be read in with slicing.
%
%   Extra information provided by ADIOSOPEN
%   ---------------------------------------
%   The ADIOS BP format stores the min/max values in each variable. 
%   The info structure therefore contains these min/max values. There is
%   practically no overhead to provide this information (along with the
%   values of all attributes) even for file sizes of several terabytes.
%
%   Please read the file COPYING in the top directory of the ADIOS source
%   distribution for information on the copyrights.
%
%   See also ADIOSOPEN, ADIOSREAD, ADIOSREADATTR, ADIOSCLOSE

%   Copyright 2009 UT-BATTELLE, LLC
%   Date: 2018/09/07
%   Author: Norbert Podhorszki <pnorbert@ornl.gov>

help adios
