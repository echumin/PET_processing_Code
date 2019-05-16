function runs = count_runs (resid)
% count_runs.m
%
% Description:
%   Calculates the number of runs in the residuals of the model fit to the data.
%
% Usage:
%   runs = count_runs (resid)
%
% Input Variables:
%   resid    The residuals of the model fit.  The sequence of the residuals 
%            should be monotonic with time.
%
% Output Variables:
%   runs     The number of sign changes in the residuals.
%
% Dependencies:
%   none
%
% Version: 0.0.3
% Created: 26 November 2007
% Modified: 20 October 2007
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


runs=0;
for k=1:length(resid)-1
  if sign(resid(k))~=sign(resid(k+1))
    runs=runs+1;
  end
end
