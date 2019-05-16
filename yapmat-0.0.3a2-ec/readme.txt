OVERVIEW

This software (which I'm currently calling "yapmat", as in "Yet Another Pet Modeling and Analysis Toolkit") is a package of functions intended for the analysis of dynamic positron emission tomography data.  The code that performs these analyses was written by me (Marc Normandin) while a graduate student in the Weldon School of Biomedical Engineering at Purdue University and working in the Department of Radiology at the Indiana University School of Medicine.  Auxiliary code (contained in the directory src/forge) written by others has been included for specific purposes, such as reading and writing images (Jimmy Shen's NIFTI utilities, http://www.rotman-baycrest.on.ca/~jimmy/NIFTI/).  All are free (as in speech and as in beer) software and may be used in accordance with their respective licenses; for the portions written by me, this would be GNU General Public License version 3 (see license.txt) or, at your option, any later version of the GPL.


INSTALLATION

This software should run properly without modification in recent versions of Matlab (http://www.mathworks.com) or GNU Octave (http://www.gnu.org/software/octave), although you may need to install the "io" package from Octave-Forge (http://octave.sf.net) for full functionality in Octave.  If you have problems or encounter an error that you suspect is a bug, please report it as such by following the directions below in the section titled BUG REPORTING.

For the time being, all of the code provided in yapmat consists of Matlab/Octave scripts, so no compiling is needed and "installation" boils down to adding the source directories to your Matlab/Octave path.  This can be accomplished by the following steps:

1. Extract the package wherever you'd like it to reside.
2. Launch Matlab (or Octave).
3. Change your working directory to the top level of the extracted package.
4. Run the script "install.m" from the Matlab/Octave prompt.

This procedure will add the package functions to your path, but it will persist only for the current session.  If you'd like to make this path addition permanent, you should perform a final step:

5. Run the command "savepath" at the Matlab/Octave prompt.

Alternatively, you can configure your startup file (startup.m or ~/.octaverc)to execute install.m each time you launch Matlab or Octave, but please note that install.m must remain in its original position or the source directories will not be found.

If you want to uninstall yapmat, just run the uninstall.m script from your Matlab/Octave prompt.  If you'd like to remove it for good, follow this up with the "savepath" command.


DOCUMENTATION

Functions are documented in the conventional Matlab/Octave fashion, that is, you can find information by typing "help function_name" at the Matlab/Octave prompt, replacing function_name with whatever code you need to learn more about.

Several demonstrations are included in the demos/ directory.  While these are not exhaustive tutorials, they showcase some of the most commonly used features.  If you're new to yapmat (I bet you are!), these demos would be a good place to start.  Take a look at the files themselves to get a feel for how to use yapmat, then run the scripts from your Matlab/Octave prompt to see what yapmat actually does.


BUG REPORTING

If you think you've found a bug, please report it by sending me an email at normandin@ieee.org.  In the email, please include a detailed description of the problem you're encountering (including any error messages that are generated), the version of Matlab or Octave that you're using, and your operating system.  The goal is to provide enough information that I can reproduce the error, which may eventually include a request for the data that's provoking the problem; but please *do* *not* include huge data sets as attachments to you email.  This will simply cripple my email quota and put me in a foul mood, neither of which will help get your problem addressed.  If I need your data, I will ask for it and arrange for a safe way to transfer it that does not burden my inbox or my sanity.


ACKNOWLEDGEMENTS

While I wrote the code implementing the PET analysis methods included in the yapmat package, in most (currently, all) cases the techniques were conceived and published by others.  Without these concepts that I've borrowed from people smarter than me, this package would just be a collection of trivial toys (perhaps it's still just that).  So a big thank you goes out to all of the researchers whose ideas I have hijacked.  References to the publications formally documenting the methods implemented here can be found in references.txt.  I've done my best to faithfully reproduce the models as they were described, though in some cases I've made slight modifications that, in my experiences, improve the overall performance.  If these modifications, or errors in my actual implementation, have sullied the good name of a technique, I apologize and would be happy to address the issue if you file a bug report.


Thanks for choosing yapmat for your PET analysis needs!

Marc D. Normandin
normandin@ieee.org
