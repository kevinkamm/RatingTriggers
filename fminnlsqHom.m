function r=fminnlsqHom(Ut0,A,t,t0,h,PD) %#codegen
%%FMINLSQ the difference of model and market data
% generator.
%   Input:
%       Ut0 (KxK array): contains the evolution system at U_{0,t0}
%       A (KxK array): contains the market generator matrices. K is 
%                              the number of ratings
%       t0 (double): contains the year of the previous rating matrix
%       t (double): contains the year of the current rating matrix
%       h (Kx1): contains the change of measure parameters at t
%       PD (Kx1 array): contains the default probabilities at t
%   Output:
%       r (K-1x1 array): contains the difference of model and market data
h=reshape(h,[],1);
h=cat(1,h,1);
H=bsxfun(@rdivide,h',h);
Ah=(A.*H);
diag=eye(size(Ah),'logical');
Ah(diag)=0;
Ah(diag)=-sum(Ah,2);
left=Ut0*expm(Ah.*(t-t0));
right=PD(1:end-1);
r=left(1:end-1,end)-right;
end