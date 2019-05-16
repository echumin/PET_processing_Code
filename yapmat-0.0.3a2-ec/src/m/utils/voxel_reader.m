function [dyn, fstump, hdr, nii] = voxel_reader (imgfile, nframes, base);

% Reads the image data from file and stores into a variable that holds the 4D
% dynamic image data.
%
% Version: 0.0.3
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


% preprocess arguments
[fpath, fname, ftype] = fileparts (imgfile);
idxlen = length (num2str (nframes));
fstump = fname(1:end-idxlen);
if ~(strcmp (base, 'hex') || strcmp (base, 'dec'))
  error ('Invalid value for ''base'', must be ''dec'' or ''hex''.');
end

% load image data frame by frame
for i = 1:nframes
  % print progress
  fprintf ('\n  frame number %d ...', i);

  % convert dec to hex if needed
  if strcmp (base, 'hex')
    idx = lower (num2str (dec2hex (i))); % tolower faster, but not in matlab
  else % base == 'dec'
    idx = num2str (i);
  end

  % pad with zeros if needed
  while length (idx) < idxlen
    idx = strcat ('0', idx);
  end

  % read in image volume for ith frame
  nii = load_nii (fullfile (fpath, strcat (fstump, idx, ftype)));
  dyn(:,:,:,i) = nii.img;

  % print progress
  fprintf (' done.');
end

fprintf ('\nCasting dynamic image data to double precision ...');
dyn = double (dyn);
fprintf (' done.');
hdr = nii.hdr;
