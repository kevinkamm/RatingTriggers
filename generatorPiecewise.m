function [At, varargout]=generatorPiecewise(A,tMarket,t)
%%GENERATORPIECEWISE evaluates the generator in a piecewise
% fashion.
%   Input:   
%       AMarket (KxKxn array): contains the generator matrices. K is 
%                              the number of ratings and n the number of 
%                              matrices
%       tMarket (1xn array): contains the years of the rating matrices
%       t (1xN array): time points   
%   Output:
%       At (KxKxN array): contains the generator in decimals
%       varargout (cell array):
%           1 -> (1xn array): contains the times in grid nearest to tMarket
At=zeros(size(A,1),size(A,2),length(t));
tk=zeros(1,length(tMarket)+1);
for k=1:1:length(tMarket)
    tk(k+1)=find(t<=tMarket(k),1,'last');
    At(:,:,tk(k)+1:1:tk(k+1))=repmat(A(:,:,k),1,1,tk(k+1)-tk(k));
end
if nargout>1
    varargout{1}=tk(2:end);
end
end