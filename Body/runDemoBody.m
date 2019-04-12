clear; close all;
%% input files
BODY_SHAPE_ID=1;
X   =sprintf('%s/Data/vert%.3d.txt',       pwd,BODY_SHAPE_ID);
H   =sprintf('%s/Model/opa-H-lout%.3d.txt',pwd,BODY_SHAPE_ID);
mu  =sprintf('%s/Model/opa-m-lout%.3d.txt',pwd,BODY_SHAPE_ID);
dld =sprintf('%s/../win/dld.exe',pwd);

%% parameters
omg ='1e-4';
gma ='1e-4';
nsamp ='500';

%% execution
cmd=sprintf('%s %s %s %s -w %s -g %s -s -y %s -h -H',dld,X,H,mu,omg,gma,nsamp)
system(cmd); 
optpath;

