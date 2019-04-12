png=0;

fp=fopen('.optpath.bin');
N =fread(fp,  1,  'int32' );
D =fread(fp,  1,  'int32' );
M =fread(fp,  1,  'int32' );
L =fread(fp,  1,  'int32' );
T =fread(fp,D*M*L,'double');
X =fread(fp, D*N, 'double');
fclose(fp);

T =reshape(T,[D,M,L]);
X =reshape(X,[D,N]);

close;
for d=1:D

  switch(d)
    case(1) 
      R=eye(3); idx=[1,2,3];
    case(2) 
      a=1/sqrt(2); R=[a,-a,0;a,a,0;0,0,1]; idx=[3,1,2];
    case(3) 
      a=1/sqrt(2); R=[a,-a,0;a,a,0;0,0,1]; idx=[1,3,2];
  end;

  Z=R(idx,idx)*X;
  for l=1:L U(:,:,l)=R(idx,idx)*T(:,:,l); end; 
  
  for l=1:L
    Y=U(:,:,l);
    R=[0,-1,0;1,0,0;0,0,1];
    W=(R*Z)';
    Y=(R*Y)';
    plot(W(:,1),W(:,2),'b.','MarkerSize',1); hold on;
    plot(Y(:,1),Y(:,2),'r.','MarkerSize',1,'MarkerFaceColor',[1,0,0]);
    %bbox=[min(X(:,a)),max(X(:,a)),min(X(:,b)),max(X(:,b))];
    bbox=1.2*[-1,1,-1,1];

    pbaspect([1 1 1]);
    axis(bbox); axis off;
  
    if(png==1)
      fn=sprintf('Work/otw-%04d-v%d.png',l,d);
      print(fn,'-dpng');
    end;
  
    pause(0.01);
    hold off;
  end
end;
