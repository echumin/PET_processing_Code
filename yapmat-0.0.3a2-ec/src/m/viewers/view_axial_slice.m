function [] = view_axial_slice (img, slicenum, varargin)

% File: view_axial_slice.m
%
% Description:
%   Display an axial slice (or what I think is an axial slice) of an image
%   volume.
%
% Usage:
%   [ ] = view_axial_slice (IMG, SLICENUM)
%   [ ] = view_axial_slice (IMG, SLICENUM, CMAP)
%
% Inputs:
%   IMG        3D image volume, such as the output BP from mrtm2().  HINT: If 
%              the image is noisy you may have to apply thresholds, i.e., use
%              as your input "IMG .* (IMG > 0 & IMG < 10)" to display only 
%              pixels in the range 0 through 10, with values outside the range 
%              set to zero.
%   SLICENUM   An integer indicating which slice you want to display.  HINT:
%              The number should be one or greater, and no larger than the
%              output of size (IMG, 3).
%   CMAP       [Optional] Colormap used in displaying slice.  If no value is
%              provided, grayscale will be used.
%
% Outputs:
%   None, just displays the slice.
%
% Dependencies:
%   None
%
% Version: 0.0.1
% Modified: 18 August 2008
% Author: Marc Normandin (normandin@ieee.org)
%


imagesc (img(:,end:-1:1,slicenum)');
axis image;
colorbar;

if nargin == 3
  cmap = varargin{1};
  try
    colormap (cmap);
  catch
    warning ('colormap undefined, i''ll default to grayscale');
    colormap (gray);
  end
else
  warning ('colormap undefined, i''ll default to grayscale');
  colormap (gray);
end
