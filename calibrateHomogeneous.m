function [AP,AQ,hFmin,errFmin,ctimeLsqnonlin,UP,UQ]=...
                                        calibrateHomogeneous(P,tMarket,PD,muQ,muP)
%%CALIBRATEHOMOGENEOUS calibrates the model with piecewise generator using
% non-linear least squares
%   Input:
%       P (KxKxn array): contains the market rating transition probabilities. 
%                        K is the number of ratings and n the number of 
%                        matrices
%       tMarket (1xn array): contains the years of the rating matrices
%       PD (Kxn array): contains the default probabilities
%   Output:
%       AP (KxKxn array): piecewise generator under P
%       AQ (KxKxn array): piecewise generator under Q
%       hFmin (Kxn array): calibrated change of measure parameters
%       errFmin (nx1 array): error of least square norm
%       ctimeLsqnonlin (nx1 array): computational time of nonlinear LSQ

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

if muP~=0
                
% lower and upper bound for x=[A,h], such that x(:,1:end-1)=A, x(:,end)=h
% we will force h(end)=1 with the bounds
lb=1e-4.*ones(size(P,1),size(P,2)+1);lb(end,end)=1;
bndDiagInd=cat(2,eye(size(P,1),'logical'),zeros(size(P,1),1,'logical'));
lb(bndDiagInd)=-10;
lb(end,1:end-1)=0;
ub=1000.*ones(size(P,1),size(P,2)+1);
if tMarket(end)>1
ub(:,end)=5;
lb(1:end-1,end)=1e-1;
else
ub(:,end)=2;
lb(1:end-1,end)=1e-2;
end
ub(end,end)=1;
ub(bndDiagInd)=0;
ub(end,1:end-1)=0;

UPt0=eye(size(P,1));
UQt0=eye(size(P,1));
AP=zeros(size(P));
AQ=zeros(size(P));
UP=zeros(size(P));
UQ=zeros(size(P));
hFmin=zeros(size(P,1),size(P,3));
errFmin=zeros(size(P,3),1);
ctimeLsqnonlin=zeros(size(P,3),1);
time=cat(2,0,tMarket);
diagInd=eye(size(P,1),'logical');
for kk=1:1:length(tMarket)
    % Find generator under P
    temp=logm(UPt0^(-1)*P(:,:,kk))./(time(kk+1)-time(kk));
    % Adjust result to be generator
    temp(diagInd)=0;
%     temp=abs(temp);
    temp(temp<0)=0;
    temp(diagInd)=-sum(temp,2);
    AP(:,:,kk)=temp;
    
    % Iterate evolution system under P
    UPt0=UPt0*expm(AP(:,:,kk).*(time(kk+1)-time(kk)));
    UP(:,:,kk)=UPt0;
    
    % Find change of measure parameters h
    x0=lb;
    x0(:,1:end-1)=AP(:,:,kk);
    ticLSQ=tic;
    [xmin,err] = fmincon(@(x)fmin(x,AP(:,:,kk),PD(:,kk),UQt0,time(kk+1),time(kk),muQ,muP),...
                 x0,...
                 [],[],...
                 [],[],...
                 lb,ub,...
                 @(x)nonlinconstr(x,UPt0,UQt0,time(kk+1),time(kk)),...
                 opts);
    fprintf('Error of optimization: %3.3g\n',err)
    temp=xmin(:,1:end-1);
    temp(diagInd)=0;
%     temp=abs(temp);
    temp(temp<0)=0;
    temp(diagInd)=-sum(temp,2);
    AP(:,:,kk)=temp;
    fprintf('Error of final adjustment of AP: %3.3g\n',...
        sum(abs(squeeze(AP(:,:,kk))-xmin(:,1:end-1))./(size(AP,1)^2),'all'))
    hmin=xmin(:,end);
    ctimeLsqnonlin(kk)=toc(ticLSQ);
    % Add h(end+1)=1 for absorbing default
    
    % Store results
    hFmin(:,kk)=hmin;
    errFmin(kk)=err;
    
    % Build generator under Q
    Ah=(AP(:,:,kk).*(hmin'./hmin));
    Ah(diagInd)=0;
    Ah(diagInd)=-sum(Ah,2);
    AQ(:,:,kk)=Ah;
    
    % Iterate evolution system under Q
    UQt0=UQt0*expm(Ah.*(time(kk+1)-time(kk)));
    UQ(:,:,kk)=UQt0;
    
    % Outputs
%     UQt0(:,end)
%     PD(:,kk)
%     err
end
else
% lower and upper bound for x=[A,h], such that x(:,1:end-1)=A, x(:,end)=h
% we will force h(end)=1 with the bounds
lb=1e-4.*ones(size(P,1)-1,1);
ub=inf.*ones(size(lb));

UPt0=eye(size(P,1));
UQt0=eye(size(P,1));
AP=zeros(size(P));
AQ=zeros(size(P));
UP=zeros(size(P));
UQ=zeros(size(P));
hFmin=zeros(size(P,1),size(P,3));
errFmin=zeros(size(P,3),1);
ctimeLsqnonlin=zeros(size(P,3),1);
time=cat(2,0,tMarket);
diagInd=eye(size(P,1),'logical');
for kk=1:1:length(tMarket)
    % Find generator under P
    temp=logm(UPt0^(-1)*P(:,:,kk))./(time(kk+1)-time(kk));
    % Adjust result to be generator
    temp(diagInd)=0;
%     temp=abs(temp);
    temp(temp<0)=0;
    temp(diagInd)=-sum(temp,2);
    AP(:,:,kk)=temp;
    
    % Iterate evolution system under P
    UPt0=UPt0*expm(AP(:,:,kk).*(time(kk+1)-time(kk)));
    UP(:,:,kk)=UPt0;
    
    % Find change of measure parameters h
    x0=lb;
    ticLSQ=tic;
    [xmin,err] = fmincon(@(x)fmin(x,AP(:,:,kk),PD(:,kk),UQt0,time(kk+1),time(kk),muQ,muP),...
                 x0,...
                 [],[],...
                 [],[],...
                 lb,ub,...
                 @(x)nonlinconstr(x,UPt0,UQt0,time(kk+1),time(kk)),...
                 opts);
    fprintf('Error of optimization: %3.3g\n',err)

    hmin=[xmin;1];
    ctimeLsqnonlin(kk)=toc(ticLSQ);
    % Add h(end+1)=1 for absorbing default
    
    % Store results
    hFmin(:,kk)=hmin;
    errFmin(kk)=err;
    
    % Build generator under Q
    Ah=(AP(:,:,kk).*(hmin'./hmin));
    Ah(diagInd)=0;
    Ah(diagInd)=-sum(Ah,2);
    AQ(:,:,kk)=Ah;
    
    % Iterate evolution system under Q
    UQt0=UQt0*expm(Ah.*(time(kk+1)-time(kk)));
    UQ(:,:,kk)=UQt0;
    
    % Outputs
%     UQt0(:,end)
%     PD(:,kk)
%     err
end
end