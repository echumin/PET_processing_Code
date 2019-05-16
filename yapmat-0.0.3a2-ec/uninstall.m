% At the Matlab/Octave prompt, change your working directory here (i.e., 
% /path/to/package/yapmat-x.y.z) and then run this script.  That will add 
% the yapmat functions to your path for this session.  To make this 
% "installation" permanent, you'll want to save your path with the savepath 
% function.
%
% For more information, see readme.txt.
fprintf ('\nWelcome to Yet Another PET Modeling and Analysis Toolkit.');
rmpath (genpath (fullfile (pwd, 'src')));
fprintf ('\n\nYou''ve successfully uninstalled yapmat (version 0.0.3).');
fprintf ('\nTo make this ''uninstallation'' permanent, enter');
fprintf ('\n''savepath'' now at your Matlab (or Octave) prompt.');
fprintf ('\n\nSo long!\n\n')
