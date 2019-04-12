clear; close all;
drc='Model';L=40; K=39; png=0;

%% read input files
ct=0;
for l=1:L
  x=load(sprintf('Data/hand%.3d.txt',l));
  if(ct==0) [M,D]=size(x); Xorig=zeros(M,D,L); end;
  Xorig(:,:,l)=x; 
  ct=ct+1;
end;

for l=1:L
  sprintf('l=%.3d\n',l)
  %% leave-one-out
  X=Xorig; X(:,:,l)=[];
  
  %% computing eigenshapes
  X=reshape(X,[M*D,L-1])'; 
  m=(sum(X)/size(X,1))';
  [B,Lmd]=eig(cov(X'));
  A=X'*B;
  G=A*diag(1./sqrt(diag(A'*A)));
  Lmd=Lmd/(L-1);
  
  dlmwrite(sprintf('%s/H-lout%.3d.txt',drc,l),G*sqrt(Lmd),'\t');
  dlmwrite(sprintf('%s/m-lout%.3d.txt',drc,l),m,          '\t');
end;

