@echo off

set CompilerFlags=/Zi /nologo /W4 /c /D_CRT_SECURE_NO_DEPRECATE /DGRAPH_WIN32 /IC:\dislin
set LinkerFlags=/incremental:no /NOLOGO /opt:ref /DEBUG  /DEFAULTLIB:gdi32.lib user32.lib
REM 5.1 is minimum for x86, 5.2 is min for x64/ARM 

if "%1" == "" (set FileName=FDex) else (set FileName=%1)

IF NOT EXIST .\build mkdir build


IF EXIST *.pdb del *.pdb
IF EXIST *.ilk del *.ilk
IF EXIST *.obj del *.obj

call "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" x64
cl %CompilerFlags%  %FileName%.c /Fobuild\ /Fdbuild\
cd build
link %FileName%    C:\dislin\disvc.lib %LinkerFlags%
COPY /Y %FileName%.exe ..\%FileName%.exe

REM pause
