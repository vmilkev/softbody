@ECHO OFF
echo Running DEM simulation study:
%~dp0\bin\dem.exe param_data\param.dat sample\sample.smpl
echo simulation study is completed
cmd /k