clear all; close all; fclose('all'); rng(0);
availGPU=gpuDeviceCount("available");
if availGPU > 0
    device='gpu';
    gpu=gpuDevice(1);
else
    device='cpu';
end
% don't change, need multiprocessing for coefficients and multithreading
% for Magnus
delete(gcp('nocreate'));parpool('threads'); % multithreading
testSet={{'Fitch','2014',[1,3,6,12]./12,1,1},...
         {'Fitch','2014',[12]./12,1,1},...
         {'SP','2020',[12]./12,1,1},...
         {'Moody','2017',[1,3,6,12]./12,1,0},...
         {'SPlt','2020',[1,3,5,7],1,6}};
% testSet=testSet(end);
testSet=testSet(1);
for ts=1:1:length(testSet)
close('all'); fclose('all');
%gpuDevice(1);
%rng('default');
currTestSetCell=testSet(ts);
currTestSet=currTestSetCell{1};
%% Choose data set
disp('Loading data')
% Rating Matrices
ratingAgency=currTestSet{1};
ratingYear=currTestSet{2};
ratingMatrixFolder='RatingMatrix';
ratingMatrixYears=currTestSet{3};
ratingMatrixDataset=currTestSet{4};
[Padjusted,A,Pmarket,ratings]=ratingMatrixLoader(ratingMatrixFolder,...
                                                 ratingAgency,...
                                                 ratingYear,...
                                                 ratingMatrixDataset,...
                                                 ratingMatrixYears);
% Default Data
defaultProbabilityFolder='DefaultProbability';
defaultProbabilityYears=ratingMatrixYears;
defaultAgency=ratingAgency;
defaultProbabilityDataset=currTestSet{5};
[PD,~,~]=defaultProbabilityLoader(defaultProbabilityFolder,...
                                  defaultAgency,...
                                  defaultProbabilityDataset,...
                                  defaultProbabilityYears);
disp('Done loading data')
%% Calibration PHCTMC
muQ=1;
muP=1;%set to zero to remove generator calibration under P
tMarket=ratingMatrixYears;
disp('Start optimization')
ticFmin=tic;
[AP,AQ,hFmin,errFmin,ctimeLsqnonlin,UPcal,UQcal]=calibrateHomogeneous(Padjusted,tMarket,PD,muQ,muP);
ctimeFmin=toc(ticFmin);
fprintf('Finished optimization after %3.3f seconds\n',...
        ctimeFmin)
%% Calibration JLT
disp('Start optimization JLT')
ticFminJLT=tic;
[APJLT,AQJLT,hFminJLT,errFminJLT,ctimeLsqnonlinJLT,UPcalJLT,UQcalJLT]=...
    calibrateJLT(Padjusted(:,:,end),tMarket(end),PD(:,end));
ctimeFminJLT=toc(ticFminJLT);
fprintf('Finished optimization of JLT after %3.3f seconds\n',...
        ctimeFminJLT)
errJLTQPD=mean(abs(squeeze(UQcalJLT(:,end))-PD(:,end)),1);
fprintf('Average error for JLT default probabilities under Q: %1.3g\n',...
        mean(errJLTQPD,'all')./length(ratings));
% sympref('FloatingPointOutput',1);
% strUQcalJLT = latex(sym(100*UQcalJLT));
% strUQcalJLT = replace(replace(strUQcalJLT, ' &', '\,\% &'),'\\','\,\% \\');
% sympref('FloatingPointOutput',0);
%% Verify results by using iterated evo sys
errAnalyticPMarket=mean(abs(UPcal-Padjusted),[1,2]);
fprintf('Error for Evo Sys Analytic of rating transitions under P: %1.3f\n',...
        sum(errAnalyticPMarket,'all'));
errAnalyticQPD=mean(abs(squeeze(UQcal(:,end,:))-PD),1);
fprintf('Error for Evo Sys Analytic of default probabilities under Q: %1.3f\n',...
        mean(errAnalyticQPD,'all'));
%% Verify results by solving forward equation
T=ratingMatrixYears(end);
% time=cat(2,0,tMarket);
dt=1/(12^3);
ticEvoSysQ=tic;
[UQ,tk]=evoSys(dt,0,T,AQ,tMarket,@generatorPiecewise);
ctimeEvoSysQ=toc(ticEvoSysQ);
errEvoSysQPD=mean(abs(squeeze(UQ(:,end,tk))-PD),1);
fprintf('Error for Evo Sys of default probabilities under Q: %1.3f\n',...
        mean(errEvoSysQPD,'all'));
UQ(:,end,end)
PD(:,end)
ticEvoSysP=tic;
[UP,~]=evoSys(dt,0,T,AP,tMarket,@generatorPiecewise);
ctimeEvoSysP=toc(ticEvoSysP);
errEvoSysPMarket=mean(abs(UP(:,:,tk)-Padjusted),[1,2]);
fprintf('Error for Evo Sys of rating transitions under P: %1.3f\n',...
        mean(errEvoSysPMarket,'all'));
UP(:,:,end)
Padjusted(:,:,end)
%% Sample ICTMC
M=10000;
XP=zeros(length(ratings),tMarket(end)/dt+1,M);
XQ=zeros(length(ratings),tMarket(end)/dt+1,M);
XPProbs=zeros(length(ratings),length(ratings),tMarket(end)/dt+1);
XQProbs=zeros(length(ratings),length(ratings),tMarket(end)/dt+1);
disp('Start simulating ICTMC under P')
ticSimP=tic;
parfor rStart=1:1:length(ratings)
    simCTMC=tic;
    XP(rStart,:,:)=ssa(AP,tMarket,rStart,M,dt);
    fprintf('Finished simulation of ICTMC starting in %s after %1.3f seconds\n',...
            ratings{rStart},toc(simCTMC));
end
ctimeSimP=toc(ticSimP);
disp('Finished simulating ICTMC under P')
disp('Start simulating ICTMC under Q')
ticSimQ=tic;
parfor rStart=1:1:length(ratings)
    simCTMC=tic;
    XQ(rStart,:,:)=ssa(AQ,tMarket,rStart,M,dt);
    fprintf('Finished simulation of ICTMC starting in %s after %1.3f seconds\n',...
            ratings{rStart},toc(simCTMC));
end
ctimeSimQ=toc(ticSimQ);
disp('Finished simulating ICTMC under Q')
disp('Calculate transition probabilities from simulations')
ticSimProb=tic;
% for rStart=1:1:length(ratings)
for rCurr=1:1:length(ratings)
    ind=XP==rCurr;
    XPProbs(:,rCurr,:)=sum(ind,3)./M;
end
% end
% for rStart=1:1:length(ratings)
for rCurr=1:1:length(ratings)
    ind=XQ==rCurr;
    XQProbs(:,rCurr,:)=sum(ind,3)./M;
end
% end
fprintf('Finished calculating probabilities of simulations after %1.3fs\n',...
        toc(ticSimProb));
%% Print errors of simulation
errSimPMarket=mean(abs(XPProbs(:,:,tk)-Padjusted),[1,2]);
fprintf('Error of simulated transition probabilities compared to market under P: %1.3f\n',...
        mean(errSimPMarket,'all'));
errSimPEvoSys=mean(abs(XPProbs(:,:,tk)-UP(:,:,tk)),[1,2]);
fprintf('Error of simulated transition probabilities compared to evolution system under P: %1.3f\n',...
        mean(errSimPEvoSys,'all'));
errSimQEvoSys=mean(abs(XQProbs(:,:,tk)-UQ(:,:,tk)),[1,2]);
fprintf('Error of simulated transition probabilities compared to evolution system under Q: %1.3f\n',...
        mean(errSimQEvoSys,'all'));
errSimPUPcal=mean(abs(XPProbs(:,:,tk)-UPcal),[1,2]);
errSimQUQcal=mean(abs(XQProbs(:,:,tk)-UQcal),[1,2]);
%% Predefault Distribution
[Q,numOfDefaultsQ]=preDefaultDistribution(XQ,length(ratings));
preDefaultPlotsQ=plotPreDefault(Q(1:end-1,1:end-1),numOfDefaultsQ-M,ratings);
%%
[P,numOfDefaultsP]=preDefaultDistribution(XP,length(ratings));
preDefaultPlotsP=plotPreDefault(P(1:end-1,1:end-1),numOfDefaultsP-M,ratings);
%% Collateral and bilateral CVA
scale=1e6; % million euros
V=scale.*portfolio(0,T,dt,25,M);
day=linspace(1,T*365,365)./365;
ti=zeros(size(day));
N=floor((T/dt)+1);
t=linspace(0,T,N);
for i=1:1:length(ti)
    ti(i)=find(t<=day(i),1,'last');
end
m=0; %minimal transfer amount
r=0; %interest rate
RI=squeeze(XQ(1,:,:));
RC=squeeze(XQ(floor(length(ratings)/2),:,:));
LGDI=0.6;
LGDC=0.6;
%% without collateralization
thresholdsIUC=scale.*1000000.*ones(length(ratings),1);
thresholdsCUC=scale.*1000000.*ones(length(ratings),1);
[CUC,CIUC,CCUC,Vadjusted]=collateral(ti,t,V,m,RI,RC,thresholdsIUC,thresholdsCUC,r,length(ratings));
[cbvaUC,cdvaUC,ccvaUC]=CBVA(ti,t,Vadjusted,RI,RC,CUC,r,length(ratings),LGDI,LGDC);
cvaPlotsUC=plotCBVA(t,ti,Vadjusted,CUC,cbvaUC,cdvaUC,ccvaUC);
collateralPlotsUC=plotCollateral(t,ti,CIUC,CCUC,CUC,Vadjusted,RI,RC,thresholdsIUC,thresholdsCUC,ratings);
fprintf('Without collateralization:\n\tCBVA=%3.3f\n\tCDVA=%3.3f\n\tCCVA=%3.3f\n',...
    cbvaUC(end),cdvaUC(end),ccvaUC(end));
%% with collateralization depending on rating triggers
thresholdsIRT=zeros(length(ratings),1);
thresholdsIRT(1:floor(length(ratings)./2))=scale.*10;
thresholdsIRT(floor(length(ratings)./2)+1:end-2)=scale.*5;
thresholdsCRT=zeros(length(ratings),1);
thresholdsCRT(1:floor(length(ratings)./2))=scale.*10;
thresholdsCRT(floor(length(ratings)./2)+1:end-2)=scale.*5;
[CRT,CIRT,CCRT,Vadjusted]=collateral(ti,t,V,m,RI,RC,thresholdsIRT,thresholdsCRT,r,length(ratings));
[cbvaRT,cdvaRT,ccvaRT]=CBVA(ti,t,Vadjusted,RI,RC,CRT,r,length(ratings),LGDI,LGDC);
cvaPlotsRT=plotCBVA(t,ti,Vadjusted,CRT,cbvaRT,cdvaRT,ccvaRT);
collateralPlotsRT=...
    plotCollateral(t,ti,CIRT,CCRT,CRT,Vadjusted,RI,RC,thresholdsIRT,thresholdsCRT,ratings);
fprintf('With rating triggers:\n\tCBVA=%3.3f\n\tCDVA=%3.3f\n\tCCVA=%3.3f\n',...
    cbvaRT(end),cdvaRT(end),ccvaRT(end));
%% with perfect collateralization
thresholdsIPC=zeros(length(ratings),1);
thresholdsCPC=zeros(length(ratings),1);
[CPC,CIPC,CCPC,Vadjusted]=collateral(ti,t,V,m,RI,RC,thresholdsIPC,thresholdsCPC,r,length(ratings));
[cbvaPC,cdvaPC,ccvaPC]=CBVA(ti,t,Vadjusted,RI,RC,CPC,r,length(ratings),LGDI,LGDC);
cvaPlotsPC=plotCBVA(t,ti,Vadjusted,CPC,cbvaPC,cdvaPC,ccvaPC);
collateralPlotsPC=...
    plotCollateral(t,ti,CIPC,CCPC,CPC,Vadjusted,RI,RC,thresholdsIPC,thresholdsCPC,ratings);
fprintf('Perfect collateralization:\n\tCBVA=%3.3f\n\tCDVA=%3.3f\n\tCCVA=%3.3f\n',...
    cbvaPC(end),cdvaPC(end),ccvaPC(end));
%% Plot Simulations
ratingPlotsP=plotRatingModel(t,XP,floor(length(ratings)/2),ratings);
ratingPlotsQ=plotRatingModel(t,XQ,floor(length(ratings)/2),ratings);
%% Print Results
fileName=sprintf('%s_%s_%d_PD_%d_muQ_%1.3f_muP_%1.3f',ratingAgency,...
                          ratingYear,...
                          ratingMatrixDataset,...
                          defaultProbabilityDataset,...
                          muQ,...
                          muP);
%% Excel
outputExcel(fileName,...
       ratingAgency,ratingYear,...
       ratingMatrixYears,ratingMatrixDataset,...
       Padjusted,Pmarket,ratings,...
       defaultProbabilityYears,defaultProbabilityDataset,...
       PD,...
       muP,muQ,...
       AP,AQ,hFmin,UPcal,UQcal,...
       UP,UQ,tk,...
       dt,...
       M,XPProbs,XQProbs,...
       errFmin,errEvoSysPMarket,errEvoSysQPD,errSimPMarket,errSimPEvoSys,errSimQEvoSys,...
       errAnalyticPMarket,errAnalyticQPD,...
       ctimeLsqnonlin,ctimeEvoSysP,ctimeEvoSysQ,ctimeSimP,ctimeSimQ,...
       preDefaultPlotsQ,...
       r,m,LGDI,LGDC,...
       thresholdsIUC,thresholdsCUC,thresholdsIRT,thresholdsCRT,thresholdsIPC,thresholdsCPC,...
       collateralPlotsUC,collateralPlotsRT,collateralPlotsPC,...
       cvaPlotsUC,cvaPlotsRT,cvaPlotsPC);
%% PDF
outputPDF(fileName,...
       ratingAgency,ratingYear,...
       ratingMatrixYears,ratingMatrixDataset,...
       Padjusted,Pmarket,ratings,...
       defaultProbabilityYears,defaultProbabilityDataset,...
       PD,...
       muP,muQ,...
       AP,AQ,hFmin,UPcal,UQcal,...
       UP,UQ,tk,...
       dt,...
       M,XPProbs,XQProbs,...
       errFmin,errEvoSysPMarket,errEvoSysQPD,errSimPMarket,errSimPEvoSys,errSimQEvoSys,...
       errAnalyticPMarket,errAnalyticQPD,errSimPUPcal,errSimQUQcal,...
       ctimeLsqnonlin,ctimeEvoSysP,ctimeEvoSysQ,ctimeSimP,ctimeSimQ,...
       preDefaultPlotsQ,preDefaultPlotsP,...
       r,m,LGDI,LGDC,...
       cbvaUC,cdvaUC,ccvaUC,cbvaRT,cdvaRT,ccvaRT,cbvaPC,cdvaPC,ccvaPC,...
       thresholdsIUC,thresholdsCUC,thresholdsIRT,thresholdsCRT,thresholdsIPC,thresholdsCPC,...
       collateralPlotsUC,collateralPlotsRT,collateralPlotsPC,...
       cvaPlotsUC,cvaPlotsRT,cvaPlotsPC,...
       ratingPlotsP,ratingPlotsQ);
disp('done')
end