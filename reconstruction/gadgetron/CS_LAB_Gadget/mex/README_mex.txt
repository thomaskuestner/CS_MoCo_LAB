if Matlab 2013b or older:
Windows:
0.) first change the LINKFLAGS in "mexopts.bat":
file is located under: %appdata%\MathWorks\MATLAB\R2013b
if it is not existing:
	- then open Matlab and run the command: mex -setup 
	- choose a compiler and finish
inside the file "mexopts.bat":
append to the line starting with "set LINKFLAGS"
ACE.lib CS_LAB.lib /LIBPATH:lib\boost /LIBPATH:lib\ace /LIBPATH:lib\CS_LAB
new line may look something like this:
set LINKFLAGS=/dll /export:%ENTRYPOINT% /LIBPATH:"%LIBLOC%" libmx.lib libmex.lib libmat.lib /MACHINE:X64 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /manifest /incremental:NO /implib:"%LIB_NAME%.x" /MAP:"%OUTDIR%%MEX_NAME%%MEX_EXT%.map" ACE.lib CS_LAB.lib /LIBPATH:lib\boost /LIBPATH:lib\ace /LIBPATH:lib\CS_LAB

Linux:
0.) locate mexopts file in home directory or /usr folder
inside the file "mexopts":
append to the line starting with "set LINKFLAGS"
ACE.lib CS_LAB.lib /LIBPATH:lib\boost /LIBPATH:lib\ace /LIBPATH:lib\CS_LAB
new line may look something like this:
set LINKFLAGS=/dll /export:%ENTRYPOINT% /LIBPATH:"%LIBLOC%" libmx.lib libmex.lib libmat.lib /MACHINE:X64 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /manifest /incremental:NO /implib:"%LIB_NAME%.x" /MAP:"%OUTDIR%%MEX_NAME%%MEX_EXT%.map" ACE.lib CS_LAB.lib /LIBPATH:lib\boost /LIBPATH:lib\ace /LIBPATH:lib\CS_LAB

if Matlab 2014a or newer:
0.) Nothing to be done ;-)



1.) compile mex-file with:  mex CS_LAB_Interface.cpp -Iinclude\gadgetron -Iinclude\cuda -Isrc\FOCUSS -Iinclude\fftw -Iinclude -Llib\boost -Llib\ace -lACE -Llib\CS_LAB -lCS_LAB



2.) configure reconstruction with XML config file: "CS_LAB.xml"



3.) start with "Interface.m" to run example
	
