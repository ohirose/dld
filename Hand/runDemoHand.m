clear; close all;
%% input files
HAND_ID=6;
X   =sprintf('%s/Data/hand%.3d.txt',   pwd,HAND_ID);
H   =sprintf('%s/Model/H-lout%.3d.txt',pwd,HAND_ID);
mu  =sprintf('%s/Model/m-lout%.3d.txt',pwd,HAND_ID);
dld =sprintf('%s/../win/dld.exe',pwd);

%% parameters
omg ='1e-4';
gma ='1e-4';

%% execution
cmd=sprintf('%s %s %s %s -w %s -g %s -s -h -H',dld,X,H,mu,omg,gma)
system(cmd); 
optpath;

