function runint = running_integral (del, data)
% running_integral.m
%
% Description:
%   Computes the running integral of a time series using a trapezoidal approach. 
%
% Usage:
%   runint = running_integral (del,data)
%
% Input Variables:
%   del      A vector containing the duration of the PET acquisition frames.
%   data     The data whose running integral is to be calculated.  If data is a 
%            matrix, the running integral of each column will be determined.
%
% Output Variables:
%   runint   The running integral of time-series data.
%
% Dependencies:
%   none
%
% NOTE: DATA MUST BE A COLUMN VECTOR OR AN ARRAY OF COLUMN VECTORS.  NEED TO 
%       DO A CHECK AND TRANSPOSE (WITH WARNING MESSAGE) IF A ROW VECTOR IS
%       PROVIDED AS INPUT.
%
% Version: 0.0.3
% Created: 21 November 2007
% Modified: 17 January 2008
% Author: Marc Normandin (normandin@ieee.org)

% Copyright (C) 2008  Marc D. Normandin
%
% This file is part of yapmat.
%
% yapmat is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% yapmat is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with yapmat.  If not, see <http://www.gnu.org/licenses/>.

%------------------------------------------------------------------------------

% Replicate del into an array for vectorization if data contains multiple TACs
if size(data,2)~=size(del,2)
  tmp=del;
  for k=1:(size(data,2)-size(del,2))/size(del,2)
    del=[del tmp];
  end
  clear tmp;
end

% Loop over frames to compute running integral (assumes zero IC)
%  for m=1:size(del,1)
%    rint(m,:)=sum(del(1:m,:).*data(1:m,:),1)-del(m,:).*data(m,:)/2;
%  end

% Compute running integral (assumes zero initial condition)
runint=cumsum(del.*data,1)-del.*data/2;
