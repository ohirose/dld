clear; close all;
M=58; D=2; N=40; T=6;
K=50;

X=zeros(M,D,T,N);
for n=1:N for t=1:T
  fnm=sprintf('Data/%.2d-%dm-opa.txt',n,t);
  fnf=sprintf('Data/%.2d-%df-opa.txt',n,t);
  if(exist(fnm)) fn=fnm; else fn=fnf; end;
  x=load(fn);
  X(:,:,t,n)=x;
end; end;
X=reshape(X,[M*D,N*T])'; 

% leave-one-person-out (leave six faces out)
for n=1:N
  del=(n-1)*T+1:n*T;
  Y=X; Y(del,:)=[]; %size(Y)
  m=(sum(Y)/size(Y,1))';
  [G,L]=eig(cov(Y));
  [L,idx]=sort(diag(L),'descend');
  G=G(:,idx);
  H=G*diag(sqrt(L)); H=H(:,1:K);
  fnH=sprintf('Model/H-lout%.3d.txt',n); dlmwrite(fnH,H,'\t');
  fnm=sprintf('Model/m-lout%.3d.txt',n); dlmwrite(fnm,m,'\t');
end;

