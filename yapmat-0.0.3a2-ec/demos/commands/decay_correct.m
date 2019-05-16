function [] = decay_correct (imgfile, delta_t, halflife, postfix, nframes, base)
%% Reads in nifti image (imgfile) and decay corrects the data by delta_t
%% minutes. Saves result to original filename appended by user-specified
%% postfix. If postfix is not provided the default is '_dkcor'. Will loop over
%% a dynamic image sequence if number of frames (nframes) and the numerical
%% base ('hex' or 'dec') are given.
%%
%% Example:
%%   >> decay_correct ('raw001.nii', 123.4, 109.8, '_corrected', 30, 'dec');
%% The above code will decay correct the sequence of 18F PET image files 
%% raw001.nii, raw002.nii, ... raw030.nii to 123.4 min earlier. The resulting
%% decay-corrected images will be saved as raw001_corrected.nii, etc.

%% check and preprocess arguments
if ~(exist ('base', 'var'))   % argument 'nframes'
  nframes = 1;
else
if ~(exist ('base', 'var'))   % argument 'base'
  if (nframes == 1)
    base = 'dec';   % meaningless assignment so refs to base don't crash
  else
    error ('Error: Argument ''base'' not specified.');
  end
else
  base = lower (base);
  switch base
    case 'dec'
      idxlen = length (num2str (nframes));
    case 'hex'
      idxlen = length (num2str (dec2hex(nframes)));
    otherwise
      error ('Invalid value for ''base'', must be ''dec'' or ''hex''.');
  end
end
if ~(exist ('imgfile', 'var'))   % argument 'imgfile'
  error ('Error: Argument ''imgfile'' not specified.');
else
  [fpath, fname, ftype] = fileparts (imgfile);
  fstump = fname(1:end-idxlen);
  idx0 = str2num (fname(end-idxlen+1:end));
  idxn = idx0 + nframes - 1;
end
if ~(exist ('delta_t', 'var'))   % argument 'delta_t'
  error ('Error: Argument ''delta_t'' not specified.');
end
if ~(exist ('postfix', 'var'))   % argument 'postfix'
  Warning ('Warning: Argument ''postfix'' not specified. Using default ''_dkcor''.');
end
if ~(exist ('halflife', 'var'))   % argument 'halflife'
  error ('Error: Argument ''halflife'' not specified.');
end

%% loop over the frames
for ii = idx0:idxn
  %% print progress
  fprintf ('Decay correcting frame number %d of %d ...', ii, idx0+nframes);
  
  %% convert dec to hex if needed
  if strcmp (base, 'hex')
    idx = lower (num2str (dec2hex (ii)));
  else % base == 'dec'
    idx = num2str (ii);
  end

  %% pad with zeros if needed
  while length (idx) < idxlen
    idx = strcat ('0', idx);
  end

  %% read in image volume for ith frame, decay correct it
  nii = load_nii (fullfile (fpath, strcat (fstump, idx, ftype)));
  img = double (nii.img);
  img .* exp(log(2) / halflife * delta_t);
  
  %% put decay-corrected image data into nifti structure, save to file
  nii.img = img;
  save_nii (nii, strcat (fstump, idx, postfix, ftype));

  %% print progress
  fprintf (' done.\n');
end