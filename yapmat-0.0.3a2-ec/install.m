% At the Matlab/Octave prompt, change your working directory here (i.e., 
% /path/to/package/yapmat-x.y.z) and then run this script.  That will add 
% the yapmat functions to your path for this session.  To make this 
% "installation" permanent, you'll want to save your path with the savepath 
% function.
%
% For more information, see readme.txt.
fprintf ('\nWelcome to Yet Another PET Modeling and Analysis Toolkit.');
addpath (genpath (fullfile (pwd, 'src')));
fprintf ('\n\nYou''ve successfully installed yapmat (version 0.0.3).');
fprintf ('\nTo make this ''installation'' permanent, enter ''savepath''');
fprintf ('\nnow at your Matlab (or Octave) prompt.');
fprintf ('\n\nIf you''re new to yapmat (and I bet you are) you''d do well to');
fprintf ('\nlook through the readme.txt file and peruse the demo(s).');
fprintf ('\n\nPlease send comments, complaints, and other flavors of');
fprintf ('\nfeedback to Marc Normandin (normandin@ieee.org).');
fprintf ('\n\nHappy analyzing!\n\n')
