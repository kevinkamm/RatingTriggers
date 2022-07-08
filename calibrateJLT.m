function [AP,AQ,hFmin,errFmin,ctimeLsqnonlin,UP,UQ]=...
                                        calibrateJLT(P,tMarket,PD)
%%CALIBRATEHOMOGENEOUS calibrates the model with piecewise generator using
% non-linear least squares
%   Input:
%       P (K x K array): contains the market rating transition probabilities. 
%                        K is the number of ratings 
%       tMarket (double): contains the year of the rating matrix
%       PD (K x 1 array): contains the default probabilities
%   Output:
%       AP (K x K  array): piecewise generator under P
%       AQ (K x K array): piecewise generator under Q
%       hFmin (K x 1 array): calibrated change of measure parameters
%       errFmin (double): error of least square norm
%       ctimeLsqnonlin (double): computational time of nonlinear LSQ

% opts = optimoptions('fmincon',...
%                     'MaxFunctionEvaluations',100000,...
%                     'MaxIterations',1000000,...
%                     'TolCon',1e-16,...
%                     'UseParallel',true,...
%                     'Display','iter');
opts = optimoptions('fmincon',...
                    'MaxFunctionEvaluations',1000000,...
                    'MaxIterations',10000000,...
                    'ConstraintTolerance',1e-6,...
                    'OptimalityTolerance',1e-10,...
                    'StepTolerance',1e-16,...
                    'UseParallel',false);
                
% lower and upper bound for h=h_1,...,h_K-1, we are not using h_K in
% calibration since it is 1
lb=1e-4.*ones(size(P,2)-1,1);
ub=100.*ones(size(P,2)-1,1);

errFmin=0;

diagInd=eye(size(P,1),'logical');

% Find generator under P
temp=logm(P)./(tMarket);
% Adjust result to be generator
temp(diagInd)=0;
%     temp=abs(temp);
temp(temp<0)=0;
temp(diagInd)=-sum(temp,2);
AP=temp;

% Iterate evolution system under P
UP=expm(AP.*(tMarket));

ticLSQ=tic;
[xmin,errFmin] = fmincon(@(x)fminJLT(x,AP,PD,tMarket),...
             (lb+ub)./2,...
             [],[],...
             [],[],...
             lb,ub,...
             [],...
             opts);
fprintf('Error of optimization: %3.3g\n',errFmin)
hFmin=[xmin;1];

ctimeLsqnonlin=toc(ticLSQ);


% Build generator under Q
AQ=diag(hFmin)*AP;
UQ=expm(AQ.*tMarket);

end