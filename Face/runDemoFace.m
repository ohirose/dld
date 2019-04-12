clear; close all;
%% input files
HUMAN_ID=22;
FACE_ID='4';
SEX='f';

X   =sprintf('%s/Data/%.2d-%s%s-opa.txt',pwd,HUMAN_ID,FACE_ID,SEX);
H   =sprintf('%s/Model/H-lout%.3d.txt',  pwd,HUMAN_ID);
mu  =sprintf('%s/Model/m-lout%.3d.txt',  pwd,HUMAN_ID);
dld =sprintf('%s/../win/dld.exe',        pwd);

%% parameters
omg ='1e-4';
gma ='1e-2';

%% execution
cmd=sprintf('%s %s %s %s -w %s -g %s -s -h -H',dld,X,H,mu,omg,gma)
system(cmd); 
optpath;

