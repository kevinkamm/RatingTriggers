function [c,ceq]=nonlinconstr(x,UPt0,UQt0,t,t0)
%%NONLINCONSTR 
%   Input:
%       x (Kx(K+1) array): x(:,1:end-1)=A is the generator under P and 
%                          x(:,end)=h is the change of measure parameter
%   Output:
%       ceq (ndarray) : equality constraints, ceq=0
%       c (ndarray): inequality constraints, c<=0

% A=x(:,1:end-1);
% h=x(:,end);
% H=bsxfun(@rdivide,h',h);
% Ah=A.*H;
% 
% diagInd=eye(size(A,1),'logical');
% Atemp=A;
% Atemp(diagInd)=0;
% Atemp(diagInd)=-sum(Atemp,2);
% Ah(diagInd)=0;
% Ah(diagInd)=-sum(Ah,2);
% 
% UPt=UPt0*expm(Atemp.*(t-t0));
% UQt=UQt0*expm(Ah.*(t-t0));
% 
% ceq=zeros(size(A,1),2);
% c=zeros((size(A,1)-1)*size(A,1),2);
c=[];
% ceq=diag(A)+sum(Atemp,2);
ceq=[];
% c=zeros(size(UQt0,2)-1,2);
% % bd=sum(UPt(:,1:end-1),1);
% % pdFactor=abs(sum(UPt(:,end))-sum(UQt(:,end)));
% c(:,1)=-sum(UQt(:,1:end-1),1)+0.5;
% c(:,2)=sum(UQt(:,1:end-1),1)-1;

% constraint for A fulfilling the rating matrix property, Jarrow, Lando and
% Turnbull 1997 Lemma 2
% UP=flip(cumsum(flip(UPt,2),2),2);
% temp=UP(1:end-1,:)-UP(2:end,:);
% c(:,1)=temp(:);

% constraint for Ah fulfilling the rating matrix property
% UQ=flip(cumsum(flip(UQt,2),2),2);
% temp=UP(1:end-1,:)-UQ(2:end,:);
% c(:,2)=temp(:);
end