function [U,varargout]=evoSys(dt,t0,T,A,tMarket,generator)%#codegen
%%EVOSYS evaluates the evolution system via explicit Euler method for given
% generator.
%   Input:
%       dt (double): time mesh size for Euler       
%       T (double): contains the finit time horizon
%       AMarket (KxKxn array): contains the market generator matrices. K is 
%                              the number of ratings and n the number of 
%                              matrices
%       tMarket (1xn array): contains the years of the rating matrices
%       h (Kxn array): contains the change of measure parameters
%       generator (function handle): has inputs AMarket,tMarket,t,h and
%                                    gives calculates [AhModel,tk]
%   Output:
%       U (KxKxN array): contains the rating probabilities in decimals
%       varargout (cell array):
%           1 -> (1xn array): contains the times in grid nearest to tMarket

N=(T-t0)/dt + 1;
t=linspace(t0,T,N);
U=zeros(size(A,1),size(A,2),length(t));
U(:,:,1)=eye(size(A,1));
[At,tk]=generator(A,tMarket,t);
for i=1:1:length(t)-1
    U(:,:,i+1)=U(:,:,i)+U(:,:,i)*At(:,:,i).*dt;
end
if nargout>1
    varargout{1}=tk;
end
end