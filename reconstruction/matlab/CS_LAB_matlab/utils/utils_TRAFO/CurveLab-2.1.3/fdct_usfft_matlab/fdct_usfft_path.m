% fdct_usfft_path:

global CURVELABPATH
global PATHNAMESEPARATOR
global PREFERIMAGEGRAPHICS
global MATLABPATHSEPARATOR
		
PREFERIMAGEGRAPHICS = 1;
Friend = computer;

if strcmp(Friend,'MAC2'),
  PATHNAMESEPARATOR = ':';
  CURVELABPATH = ['Macintosh HD:Build 802:BMIALab', PATHNAMESEPARATOR];
  MATLABPATHSEPARATOR = ';';
elseif isunix,
  PATHNAMESEPARATOR = '/';
  CURVELABPATH = [pwd, PATHNAMESEPARATOR];
  MATLABPATHSEPARATOR = ':';
elseif strcmp(Friend(1:2),'PC');
  PATHNAMESEPARATOR = '\';	  
  CURVELABPATH = [pwd, PATHNAMESEPARATOR];  
  MATLABPATHSEPARATOR = ';';
end

post = PATHNAMESEPARATOR;
p = path;
pref = [MATLABPATHSEPARATOR CURVELABPATH];
p = [p pref];

p = [p pref 'CurveCoeff' post];
p = [p pref 'USFFT' post];
p = [p pref 'Utilities' post];
p = [p pref 'Windows' post 'Meyer' post ];
p = [p pref 'Windows' post 'IteratedSine' post];

path(p);

clear p pref post
clear BMIALABPATH MATLABVERSION PATHNAMESEPARATOR
clear Friend PREFERIMAGEGRAPHICS MATLABPATHSEPARATOR
