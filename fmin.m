function r=fmin(x,Amarket,PD,UQt0,t,t0,muQ,muP)
if muP==0
    A=Amarket;
    h=[x;1];
    diagInd=eye(size(A,1),'logical');
    % A(diagInd)=0;
    % A(diagInd)=-sum(A,2);
    H=bsxfun(@rdivide,h',h);
    Ah=(A.*H);
    Ah(diagInd)=0;
    Ah(diagInd)=-sum(Ah,2);
    leftQ=UQt0*expm(Ah.*(t-t0));
    rightQ=PD;
    r=muQ*sum(abs(leftQ(1:end-1,end)-rightQ(1:end-1)).^2);
else
    A=x(:,1:end-1);
    h=x(:,end);
    diagInd=eye(size(A,1),'logical');
    % A(diagInd)=0;
    % A(diagInd)=-sum(A,2);
    H=bsxfun(@rdivide,h',h);
    Ah=(A.*H);
    Ah(diagInd)=0;
    Ah(diagInd)=-sum(Ah,2);
    leftQ=UQt0*expm(Ah.*(t-t0));
    rightQ=PD;
    r=muQ*sum(abs(leftQ(1:end-1,end)-rightQ(1:end-1)).^2)+muP*sum(abs(A-Amarket).^2,'all');
end
end