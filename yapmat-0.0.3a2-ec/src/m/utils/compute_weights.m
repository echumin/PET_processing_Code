function [w1] = compute_weights (PET, ts, del, lambda, q)%,nvox)
%function [w1,w2] = compute_weights (PET, ts, del, lambda, q)%,nvox)
% function [w,sd]=computeWeights(PET,ts,del,tracer,q)
% function w=computeWeights(PET,nPETvx,REF,nREFvx,ts,del,tracer,q);    % account for region sizes, uncertainty in input
% computeWeights.m
%
% Description:
%   Calculates the weights to be applied during parameter estimation.
%
% Usage:
%   [wPET,wRintPET]=computeWeights(PET,frames,model,tracer)
%
% Input Variables:
%   PET      Time activity curve from the region of interest.
%   ts       A vector holding the start times of the PET acquisition frames.
%   del      A vector containing the duration of the PET acquisition frames.
%   model    Noise model to be used to calculate the weights.  All model 
%            formulations set weight proportional to inverse of the variance.
%              (1) sigma(i)=q*sqrt(PET(i)/delt(i))
%              (2) sigma(i)=q1+q2*sqrt(PET(i)/delt(i))
%              (3) sigma(i)=q1+q2*(PET(i)/delt(i))^q3
%            Model (2) is a degenerate of (3) with q3=0.5.  Similarly, (1) is a
%            degenerate of (2) with q1=0;  Model (1) is appropriate if assuming
%            that number of detected events (counts) is a Poisson random 
%            variable.  
%   tracer   Indicates the tracer that was used in the PET study.  The tracer 
%            (or more specifically, the tracer isotope) is needed to decay the 
%            PET data for proper calculation of sample variance, and thus 
%            weighting.
%
% Output Variable:
%   wPET     Vector with square roots of the weights to be used in parameter
%            estimation for the time series data.
%   wRintPET Vector with square roots of the weights to be used in parameter
%            estimation for the running integral of the time series data.
%
% Dependencies: 
%   running_integral.m
%
% Version: 0.0.3
% Created: 21 November 2007
% Modified: 20 October 2008
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

% Calculate frame midpoints, in minutes
midpt=(ts+del/2);

% Determine the tracer isotope decay constant
if ischar (lambda)
  switch lambda
    case '11C'
      lambda=log(2)/20.3;
    case '18F'
      lambda=log(2)/109.8;
    otherwise
      warning('Isotope not supported.  I''ll perform unweighted fitting.');
      w=ones(length(ts),1);
      return;
  end
end

% Default noise model is (1), counts are Poisson-like RV
%q=[0 1 0.5]';

% Calculate weights for PET data
dkPET=PET.*exp(-lambda.*midpt);
sd1=abs(q(1)+q(2)*(dkPET./del).^q(3));   % abs() in case "negative" activity
if sum(sd1==0)~=0
  sd1(find(sd1==0))=1e-6;
end
w1=1./sd1;

% Calculate weights for running integral of PET data
%counts=running_integral(del,dkPET);  % *nPETvx;
%delt=cumsum(del);
%sd2=abs(q(1)+q(2)*(counts./delt).^q(3));  % abs() in case of "negative" counts
%w2=1./sd2;



% NEED TO ACCOUNT FOR UNCERTAINTY IN INPUT FUNCTION
% dkREF=PET.*exp(-lambda.*midpt);
% REFcts=REF.*exp(-lambda.*midpt)*nREFvx;
% sdREF=q(1)+q(2)*(REFcts./delt).^q3
% sdTOT = (sdPET+sdREF) OR sqrt((sdPET^2)+(sdREF^2)) ?????
