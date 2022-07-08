function outputPDF(fileName,...
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
               shortRate,m,LGDI,LGDC,...
               cbvaUC,cdvaUC,ccvaUC,cbvaRT,cdvaRT,ccvaRT,cbvaPC,cdvaPC,ccvaPC,...
               thresholdsIUC,thresholdsCUC,thresholdsIRT,thresholdsCRT,thresholdsIPC,thresholdsCPC,...
               collateralPlotsUC,collateralPlotsRT,collateralPlotsPC,...
               cvaPlotsUC,cvaPlotsRT,cvaPlotsPC,...
               ratingPlotsP,ratingPlotsQ)
fclose('all');
numberDict={'One','Two','Three','Four','Five','Six','Seven','Eight','Nine'};
identifierRating=[ratingAgency(1),...
                  char(65+str2num(ratingYear(end-1:end))),...
                  char(65+ratingMatrixDataset)];
identifierDefault=char(65+defaultProbabilityDataset);
identifier=[identifierRating,identifierDefault];

if length(ratings)>7
    landscape=1;
else
    landscape=0; 
end

picType='eps';
saveParam='epsc';

root=[pwd, '\' ,'Results'];
pdfRoot=[root,'\','Pdf'];
tempPath=[pdfRoot,'\','temp'];
copyPath=[pdfRoot,'\',fileName];
templatePath=[tempPath,'\','template', '.','tex'];
outputFilePath=[copyPath,'\','template','.','pdf';...
            copyPath,'\','template','.','tex'];
copyFilePath=[copyPath,'\',fileName,'.','pdf';...
              copyPath,'\',fileName,'.','tex'];
inputPath=[tempPath,'\','input','.','tex'];
          
mkDir(pdfRoot);
mkDir(tempPath);
cleanDir(tempPath,{'template.tex'});
delDir(copyPath)
mkDir(copyPath)
delFile(inputPath);

inputFile=fopen(inputPath,'a');
%% Head
fprintf(inputFile,...
        '\\section{ICTMC}\n');
%% Calibration Parameters
fprintf(inputFile,...
        '\\subsection{Calibration}\n');
% Input Param
[latexFilePath,latexCommand]=saveCalibrationParam('Calibration');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:1:length(latexCommand)
    fprintf(inputFile,...
        '%s%%\n',latexCommand{iLC});
end
% Output Generators
[latexFilePath,latexCommand]=saveCalibration('Calibration');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:1:length(latexCommand)
if landscape && iLC>1
    fprintf(inputFile,...
        '\\begin{landscape}\n');
end
fprintf(inputFile,...
        '\\begin{table}[h!]\n');
if landscape && iLC>1
    fprintf(inputFile,...
        '\\resizebox{\\linewidth}{!}{%%\n');
end
fprintf(inputFile,...
        '\t%s%%\n',latexCommand{iLC});
if landscape && iLC>1
    fprintf(inputFile,...
        '}');
end
fprintf(inputFile,...
        '\\centering\n');
fprintf(inputFile,...
        '%sCaption\n',latexCommand{iLC});
fprintf(inputFile,...
        '\\label{tab:calibration%i}\n',iLC);
fprintf(inputFile,...
        '\\end{table}\n');
if landscape && iLC>1
    fprintf(inputFile,...
        '\\end{landscape}\n');
end
end
%% Errors
fprintf(inputFile,...
        '\\section{Errors}\n');
[latexFilePath,latexCommand]=saveErrors('Errors');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:1:length(latexCommand)
% if landscape
%     fprintf(inputFile,...
%         '\\begin{landscape}\n');
% end
fprintf(inputFile,...
        '\\begin{table}[h!]\n');
fprintf(inputFile,...
        '\t%s\n',latexCommand{iLC});
fprintf(inputFile,...
        '\\centering\n');
fprintf(inputFile,...
        '%sCaption\n',latexCommand{iLC});
fprintf(inputFile,...
        '\\label{tab:errors%i}\n',iLC);
fprintf(inputFile,...
        '\\end{table}\n');
% if landscape
%     fprintf(inputFile,...
%         '\\end{landscape}\n');
% end
end
%% Computational times
fprintf(inputFile,...
        '\\subsection{Computational Times}\n');
[latexFilePath,latexCommand]=saveCompTimes('CompTimes');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:1:size(latexCommand,1)
    fprintf(inputFile,...
            '\t%s\\hfill\\\\\n',latexCommand{iLC});
end
%% Figures
fprintf(inputFile,...
        '\\section{Figures}\n');
%% Plots of Predefault
fprintf(inputFile,...
        '\\subsection{Pre-default under $\\mathbb{Q}$}\n');
latexFilePath=saveFigures(preDefaultPlotsQ,'PrePDQ',['PrePDQ','_',identifier]);
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
    fprintf(inputFile,...
        '\\subsection{Pre-default under $\\mathbb{P}$}\n');
latexFilePath=saveFigures(preDefaultPlotsP,'PrePDP',['PrePDP','_',identifier]);
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
%% Plots of collateral
fprintf(inputFile,...
        '\\subsection{Collateral and CXVA, $X=B,D,C$}\n');
fprintf(inputFile,...
        '\\subsubsection{Un-collateralized}\n');
[latexFilePath,latexCommand]=saveCollateralParam(...
                                thresholdsIUC,thresholdsCUC,cbvaUC(end),cdvaUC(end),ccvaUC(end),'UC','CollateralUC');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:1:length(latexCommand)-2
    fprintf(inputFile,...
            '\t%s\\hfill\\\\\n',latexCommand{iLC});
end
for iLC=length(latexCommand)-1:1:length(latexCommand)
if landscape
    fprintf(inputFile,...
        '\\begin{landscape}\n');
end
fprintf(inputFile,...
        '\\begin{table}[h!]\n');
if landscape
    fprintf(inputFile,...
        '\\resizebox{\\linewidth}{!}{%%\n');
end
fprintf(inputFile,...
        '\t%s%%\n',latexCommand{iLC});
if landscape
    fprintf(inputFile,...
        '}');
end
fprintf(inputFile,...
        '\\centering\n');
fprintf(inputFile,...
        '%sCaption\n',latexCommand{iLC});
fprintf(inputFile,...
        '\\label{tab:thresholds%i}\n',iLC);
fprintf(inputFile,...
        '\\end{table}\n');
if landscape
    fprintf(inputFile,...
        '\\end{landscape}\n');
end
end
latexFilePath=saveFigures(collateralPlotsUC,'CollateralUC',['CollateralUC','_',identifier]);
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
% Plots of cva
fprintf(inputFile,...
        '\\paragraph*{CVA}\n');
latexFilePath=saveFigures(cvaPlotsUC,'CVAUC',['CVAUC','_',identifier]);
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
fprintf(inputFile,...
        '\\subsubsection{Collateralized with rating triggers}\n');
[latexFilePath,latexCommand]=saveCollateralParam(...
                                thresholdsIRT,thresholdsCRT,cbvaRT(end),cdvaRT(end),ccvaRT(end),'RT','CollateralRT');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:1:length(latexCommand)-2
    fprintf(inputFile,...
            '\t%s\\hfill\\\\\n',latexCommand{iLC});
end
for iLC=length(latexCommand)-1:1:length(latexCommand)
if landscape
    fprintf(inputFile,...
        '\\begin{landscape}\n');
end
fprintf(inputFile,...
        '\\begin{table}[h!]\n');
if landscape
    fprintf(inputFile,...
        '\\resizebox{\\linewidth}{!}{%%\n');
end
fprintf(inputFile,...
        '\t%s%%\n',latexCommand{iLC});
if landscape
    fprintf(inputFile,...
        '}');
end
fprintf(inputFile,...
        '\\centering\n');
fprintf(inputFile,...
        '%sCaption\n',latexCommand{iLC});
fprintf(inputFile,...
        '\\label{tab:thresholds%i}\n',iLC);
fprintf(inputFile,...
        '\\end{table}\n');
if landscape
    fprintf(inputFile,...
        '\\end{landscape}\n');
end
end
latexFilePath=saveFigures(collateralPlotsRT,'CollateralRT',['CollateralRT','_',identifier]);
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
% Plots of cva
fprintf(inputFile,...
        '\\paragraph*{CVA}\n');
latexFilePath=saveFigures(cvaPlotsRT,'CVART',['CVART','_',identifier]);
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
fprintf(inputFile,...
        '\\subsubsection{Perfectly collateralized}\n');
[latexFilePath,latexCommand]=saveCollateralParam(...
                                thresholdsIPC,thresholdsCPC,cbvaPC(end),cdvaPC(end),ccvaPC(end),'PC','CollateralPC');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:1:length(latexCommand)-2
    fprintf(inputFile,...
            '\t%s\\hfill\\\\\n',latexCommand{iLC});
end
for iLC=length(latexCommand)-1:1:length(latexCommand)
if landscape
    fprintf(inputFile,...
        '\\begin{landscape}\n');
end
fprintf(inputFile,...
        '\\begin{table}[h!]\n');
if landscape
    fprintf(inputFile,...
        '\\resizebox{\\linewidth}{!}{%%\n');
end
fprintf(inputFile,...
        '\t%s%%\n',latexCommand{iLC});
if landscape
    fprintf(inputFile,...
        '}');
end
fprintf(inputFile,...
        '\\centering\n');
fprintf(inputFile,...
        '%sCaption\n',latexCommand{iLC});
fprintf(inputFile,...
        '\\label{tab:thresholds%i}\n',iLC);
fprintf(inputFile,...
        '\\end{table}\n');
if landscape
    fprintf(inputFile,...
        '\\end{landscape}\n');
end
end
latexFilePath=saveFigures(collateralPlotsPC,'CollateralPC',['CollateralPC','_',identifier]);
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
% Plots of cva
fprintf(inputFile,...
        '\\paragraph*{CVA}\n');
latexFilePath=saveFigures(cvaPlotsPC,'CVAPC',['CVAPC','_',identifier]);
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
%% Plots of Rating model trajectories
fprintf(inputFile,...
        '\\subsection{Rating Model trajectories}\n');
fprintf(inputFile,...
        '\\subsubsection{Under $\\mathbb{P}$}\n');
latexFilePath=saveFigures(ratingPlotsP,'RatingPlotsP',['RPP','_',identifier]);
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
fprintf(inputFile,...
        '\\subsubsection{Under $\\mathbb{Q}$}\n');
latexFilePath=saveFigures(ratingPlotsQ,'RatingPlotsQ',['RPQ','_',identifier]);
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
%% Appendix
fprintf(inputFile,...
        '\\appendix\n');
% Rating Matrices
fprintf(inputFile,...
        '\\section{Rating Matrices}\n');
[latexFilePath,latexCommand]=saveRatingMatrices('RatingMatrices');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:1:length(latexCommand)
if landscape
    fprintf(inputFile,...
        '\\begin{landscape}\n');
end
fprintf(inputFile,...
        '\\begin{table}[h!]\n');
if landscape
    fprintf(inputFile,...
        '\\resizebox{\\linewidth}{!}{%%\n');
end
fprintf(inputFile,...
        '\t%s%%\n',latexCommand{iLC});
if landscape
    fprintf(inputFile,...
        '}');
end
fprintf(inputFile,...
        '\\centering\n');
fprintf(inputFile,...
        '%sCaption\n',latexCommand{iLC});
fprintf(inputFile,...
        '\\label{tab:ratingMatrix%i}\n',iLC);
fprintf(inputFile,...
        '\\end{table}\n');
if landscape
    fprintf(inputFile,...
        '\\end{landscape}\n');
end
end
% Default Probabilities
fprintf(inputFile,...
        '\\section{Default Probabilities}\n');
[latexFilePath,latexCommand]=saveDefaultProbabilities('DefaultProbabilities');
fprintf(inputFile,...
        '\t\\input{%s}\n',changeSlash(latexFilePath));
for iLC=1:1:length(latexCommand)
fprintf(inputFile,...
        '\\begin{table}[h!]\n');
fprintf(inputFile,...
        '\t%s\n',latexCommand{iLC});
fprintf(inputFile,...
        '\\centering\n');
fprintf(inputFile,...
        '%sCaption\n',latexCommand{iLC});
fprintf(inputFile,...
        '\\label{tab:defaultProbabilities%i}\n',iLC);
fprintf(inputFile,...
        '\\end{table}\n');
end
fclose(inputFile);
% Compile Latex
currFolder=cd(tempPath);
str1=sprintf('pdflatex %s',templatePath);
[out,~]=system(str1);
cd(currFolder);

copyfile(tempPath,copyPath);
% Renaming files
for iFile=1:1:size(copyFilePath,1)
    movefile(outputFilePath(iFile,:),copyFilePath(iFile,:))
end
% Delete auxiliary latex files
delete([copyPath,'\','template*.*']);
%% Latex generating functions
%% Collateral
    function [latexFilePath,latexCommand]=saveCollateralParam(...
             thresholdI,thresholdC,cbva,cdva,ccva,colCase,saveAt)
        if strcmp(saveAt, '')
            latexFilePath='collateralParam.tex';
        else
            latexFilePath=[saveAt,'\','collateralParam.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        latexCommand={['\shortRate',colCase,identifier];...
                      ['\LGDI',colCase,identifier];
                      ['\LGDC',colCase,identifier];
                      ['\simulations',colCase,identifier];
                      ['\cbva',colCase,identifier];
                      ['\cdva',colCase,identifier];
                      ['\ccva',colCase,identifier];
                      ['\tresholdsI',colCase,identifier];
                      ['\tresholdsC',colCase,identifier]};
        file = fopen([tempPath,'\',latexFilePath],'a'); 
        fprintf(file,'\\newcommand{%s}{%%\n',latexCommand{1});
        fprintf(file,'$r=%1.3g$%%\n',shortRate);
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%svalue}{%%\n',latexCommand{1});
        fprintf(file,'%1.3g%%\n',shortRate);
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%s}{%%\n',latexCommand{2});
        fprintf(file,'$\\mathrm{LGD}_I=%1.3g$%%\n',LGDI);
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%svalue}{%%\n',latexCommand{2});
        fprintf(file,'%1.3g%%\n',LGDI);
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%s}{%%\n',latexCommand{3});
        fprintf(file,'$\\mathrm{LGD}_C=%1.3g$%%\n',LGDC);
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%svalue}{%%\n',latexCommand{3});
        fprintf(file,'%1.3g%%\n',LGDC);
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%s}{%%\n',latexCommand{4});
        fprintf(file,'$M=%d$%%\n',M);
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%svalue}{%%\n',latexCommand{4});
        fprintf(file,'%d%%\n',M);
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%s}{%%\n',latexCommand{5});
        fprintf(file,'$\\mathrm{CBVA}=%3.3g$%%\n',cbva);
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%svalue}{%%\n',latexCommand{5});
        fprintf(file,'%3.3g%%\n',cbva);
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%s}{%%\n',latexCommand{6});
        fprintf(file,'$\\mathrm{CDVA}=%3.3g$%%\n',cdva);
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%svalue}{%%\n',latexCommand{6});
        fprintf(file,'%3.3g%%\n',cdva);
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%s}{%%\n',latexCommand{7});
        fprintf(file,'$\\mathrm{CCVA}=%3.3g$%%\n',ccva);
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%svalue}{%%\n',latexCommand{7});
        fprintf(file,'%3.3g%%\n',ccva);
        fprintf(file,'}%%\n');
        header='';
        for i=1:1:length(ratings)
            if i<length(ratings)
                header=[header,sprintf('%s & ',ratings{i})];
            else
                header=[header,sprintf('%s ',ratings{i})];
            end
        end
        
        fprintf(file,'\\newcommand{%sHeader}{%%\n',latexCommand{8});
        fprintf(file,header);
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%sBody}{%%\n',latexCommand{8});
        fprintf(file,mat2TableBody({},thresholdI','g',2,2,''));
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%s}{%%\n',latexCommand{8});
        fprintf(file,'\\begin{tabular}{|*{%d}{c}|}\n',length(ratings));
        fprintf(file,'\\hline\n');
        fprintf(file,[header,'\\\\']);
        fprintf(file,mat2TableBody({},thresholdI','g',2,2,''));
        fprintf(file,'\\\\\\hline\n');
        fprintf(file,'\\end{tabular}\n');
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%sCaption}{%%\n',latexCommand{8});
        fprintf(file,'\\caption{Investor''s thresholds}');
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%sHeader}{%%\n',latexCommand{9});
        fprintf(file,header);
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%sBody}{%%\n',latexCommand{9});
        fprintf(file,mat2TableBody({},thresholdC','g',2,2,''));
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%s}{%%\n',latexCommand{9});
        fprintf(file,'\\begin{tabular}{|*{%d}{c}|}\n',length(ratings));
        fprintf(file,'\\hline\n');
        fprintf(file,[header,'\\\\']);
        fprintf(file,mat2TableBody({},thresholdC','g',2,2,''));
        fprintf(file,'\\\\\\hline\n');
        fprintf(file,'\\end{tabular}\n');
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%sCaption}{%%\n',latexCommand{9});
        fprintf(file,'\\caption{Counterparty''s thresholds}');
        fprintf(file,'}%%\n');
        fclose(file);
    end
%% Calibration
    function [latexFilePath,latexCommand]=saveCalibrationParam(saveAt)
        if strcmp(saveAt, '')
            latexFilePath='calibrationParam.tex';
        else
            latexFilePath=[saveAt,'\','calibrationParam.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        latexCommand={['\muP',identifier];...
                      ['\muQ',identifier]};
        file = fopen([tempPath,'\',latexFilePath],'a'); 
        fprintf(file,'\\newcommand{%s}{%%\n',latexCommand{1});
        fprintf(file,'$\\mu_{\\mathbb{P}}=%1.3g$%%\n',muP);
        fprintf(file,'}%%\n');
        fprintf(file,'\\newcommand{%s}{%%\n',latexCommand{2});
        fprintf(file,'$\\mu_{\\mathbb{Q}}=%1.3g$%%\n',muQ);
        fprintf(file,'}%%\n');
        fclose(file);
    end
    function [latexFilePath,latexCommand]=saveCalibration(saveAt)
        if strcmp(saveAt, '')
            latexFilePath='calibration.tex';
        else
            latexFilePath=[saveAt,'\','calibration.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        latexCommandTemp={['\calibrationH',identifier];...
                      ['\generatorAP',identifier];
                      ['\generatorAQ',identifier]};
        latexCommand={};
        description=ratings;
        generatorMatrices=zeros(size(AP,1),size(AP,2),size(AP,3),2);
        generatorMatrices(:,:,:,1)=AP;
        generatorMatrices(:,:,:,2)=AQ;
        file = fopen([tempPath,'\',latexFilePath],'a'); 
        for r=1:1:length(latexCommandTemp)
            for k=1:1:size(Padjusted,3)
                % Header
                if r==1
                    if k==1
                        latexCommand{end+1}=...
                        sprintf('%s',latexCommandTemp{r});
                        fprintf(file,...
                                '\\newcommand{%s}{\n',latexCommandTemp{r});
                        if size(Padjusted,3)==1
                        fprintf(file,...
                            '\\begin{tabular}{|c|*{%d}{c}|}\n',1);
                        else
                        fprintf(file,...
                            '\\begin{tabular}{|c|*{%d}{c}}\n',1);
                        end
                        fprintf(file,'\\hline\n');
                        fprintf(file,'Rating & $h_{%1.3f}$\\\\',ratingMatrixYears(k));
                        fprintf(file,'\\hline\n');
                        fprintf(file,...
                        mat2TableBody(description,...
                                      hFmin(:,k),...
                                      'g',2,2,''));
                    else
                    if k==size(Padjusted,3)
                        fprintf(file,...
                            '\\begin{tabular}{*{%d}{c}|}\n',1);
                    else
                        fprintf(file,...
                            '\\begin{tabular}{*{%d}{c}}\n',1);
                    end
                    fprintf(file,'\\hline\n');
                    fprintf(file,'$h_{%1.3f}$\\\\',ratingMatrixYears(k));
                    fprintf(file,'\\hline\n');
                    fprintf(file,...
                    mat2TableBody({},...
                                  hFmin(:,k),...
                                  'g',2,2,''));
                    end
                else
                    latexCommand{end+1}=...
                    sprintf('%s%s',latexCommandTemp{r},numberDict{k});
                    fprintf(file,...
                            '\\newcommand{%s%s}{\n',latexCommandTemp{r},numberDict{k});
                    fprintf(file,...
                        '\\begin{tabular}{|c|*{%d}{c}|}\n',length(ratings));
                    fprintf(file,'\\hline\n');
                    fprintf(file,'\\diagbox{From}{To} & ');
                    for i=1:1:length(ratings)
                        if i<length(ratings)
                            fprintf(file,'%s & ',ratings{i});
                        else
                            fprintf(file,'%s\\\\ ',ratings{i});
                        end
                    end
                    fprintf(file,'\\hline\n');
                    fprintf(file,...
                    mat2TableBody(description,...
                                  generatorMatrices(:,:,k,r-1),...
                                  'g',1,2,''));
                end
                if r==1
                    fprintf(file,'\\\\\\hline\n');
                    fprintf(file,...
                            '\\end{tabular}%%\n');
                    if k==size(Padjusted,3)
                        fprintf(file,...
                        '}\n');
                        fprintf(file,...
                                '\\newcommand{%sCaption}{\n',latexCommandTemp{r});
                        fprintf(file,...
                            '\\caption{Change of measure parameters}');
                        fprintf(file,'}\n');
                    end
                else
                fprintf(file,'\\\\\\hline\n');
                fprintf(file,...
                        '\\end{tabular}\n');
                fprintf(file,...
                        '}\n');
                fprintf(file,...
                        '\\newcommand{%s%sCaption}{\n',latexCommandTemp{r},numberDict{k});
                fprintf(file,...
                    '\\caption{%s at %1.3f}',...
                    replace(...
                        replace(...
                            replace(...
                                latexCommandTemp{r},...
                            '\generator',''),...
                        identifier,''),...
                    '\calibration',''),...
                    ratingMatrixYears(k));
                fprintf(file,'}\n');
                end
            end
        end
        fclose(file);
    end
%% Errors
    function [latexFilePath,latexCommand]=saveErrors(saveAt)
        if strcmp(saveAt, '')
            latexFilePath='errors.tex';
        else
            latexFilePath=[saveAt,'\','errors.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
%         latexCommandTemp={['\errorsFmin',identifier];...
%                           ['\errorsEvoSysPMarket',identifier];
%                           ['\errorsEvoSysQPD',identifier];
%                           ['\errorsSimPMarket',identifier];
%                           ['\errorsSimPEvoSys',identifier];
%                           ['\errorserrSimQEvoSys',identifier]};
        latexCommand={['\errors',identifier]};
        description={'Fmin',...
                     '$\\mathrm{U}^{\\mathbb{P}}_{\\mathrm{kol}}-\\mathrm{R}$',...
                     '$\\mathrm{U}^{\\mathbb{Q}}_{\\mathrm{kol}}e_K-\\mathrm{PD}$',...
                     '$\\mathrm{U}^{\\mathbb{P}}_{\\mathrm{ana}}-\\mathrm{R}$',...
                     '$\\mathrm{U}^{\\mathbb{Q}}_{\\mathrm{ana}}e_K-\\mathrm{PD}$',...
                     '$\\mathrm{XP}^{\\mathbb{P}}-\\mathrm{R}$',...
                     '$\\mathrm{XP}^{\\mathbb{P}}-\\mathrm{U}^{\\mathbb{P}}$',...
                     '$\\mathrm{XQ}^{\\mathbb{Q}}-\\mathrm{U}^{\\mathbb{Q}}$'};
        errors=zeros(size(Padjusted,3),6);
        errors(:,1)=errFmin;
        errors(:,2)=errEvoSysPMarket;
        errors(:,3)=errEvoSysQPD;
        errors(:,4)=errAnalyticPMarket;
        errors(:,5)=errAnalyticQPD;
        errors(:,6)=errSimPMarket;
        errors(:,7)=errSimPEvoSys;
        errors(:,8)=errSimQEvoSys;
        file = fopen([tempPath,'\',latexFilePath],'a');
        for i=1:1:size(errFmin,1)
            fprintf(file,...
                    '\\newcommand{%s}{\n',['\errorsFmin',identifier,numberDict{i},'value']);
            fprintf(file,'%3.3g',errFmin(i));
            fprintf(file,...
                    '}\n');
            fprintf(file,...
                    '\\newcommand{%s}{\n',['\errorsAnalyticPMarket',identifier,numberDict{i},'value']);
            fprintf(file,'%3.3g',errAnalyticPMarket(i));
            fprintf(file,...
                    '}\n');
            fprintf(file,...
                    '\\newcommand{%s}{\n',['\errorsAnalyticQPD',identifier,numberDict{i},'value']);
            fprintf(file,'%3.3g',errAnalyticQPD(i));
            fprintf(file,...
                    '}\n');
            fprintf(file,...
                    '\\newcommand{%s}{\n',['\errorsAnalyticPSimP',identifier,numberDict{i},'value']);
            fprintf(file,'%3.3g',errSimPUPcal(i));
            fprintf(file,...
                    '}\n');
            fprintf(file,...
                    '\\newcommand{%s}{\n',['\errorsAnalyticQSimQ',identifier,numberDict{i},'value']);
            fprintf(file,'%3.3g',errSimQUQcal(i));
            fprintf(file,...
                    '}\n');
        end
        fprintf(file,...
                '\\newcommand{%s}{\n',latexCommand{1});
        fprintf(file,...
                '\\begin{tabular}{|c|*{%d}{c}|}\n',length(ratingMatrixYears));
        fprintf(file,'\\hline\n');
        % Header
        fprintf(file,'\\diagbox{Error}{Time} & ');
        for i=1:1:length(ratingMatrixYears)
            if i<length(ratingMatrixYears)
                fprintf(file,'$%1.3f$ & ',ratingMatrixYears(i));
            else
                fprintf(file,'$%1.3f$ \\\\\n',ratingMatrixYears(i));
            end
        end
        fprintf(file,'\\hline\n');

        % Body
        fprintf(file,...
            mat2TableBody(description,...
                          errors',...
                          'g',3,3,''));
        fprintf(file,'\\\\\\hline\n');
        fprintf(file,...
                '\\end{tabular}\n');
        fprintf(file,...
                '}\n');
        fprintf(file,...
                '\\newcommand{%sCaption}{\n',latexCommand{1});
        fprintf(file,...
            '\\caption{Errors}');
        fprintf(file,...
                '}\n');
        fclose(file);
    end
%% Computational times
    function [latexFilePath,latexCommand]=saveCompTimes(saveAt)
        latexCommand={['\ctime',identifier]};
        if strcmp(saveAt, '')
            latexFilePath='ctime.tex';
        else
            latexFilePath=[saveAt,'\','ctime.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        data=[[reshape(ctimeLsqnonlin,1,[]),sum(ctimeLsqnonlin,'all')];...
              [zeros(1,length(ctimeLsqnonlin)),ctimeEvoSysP];
              [zeros(1,length(ctimeLsqnonlin)),ctimeEvoSysQ];
              [zeros(1,length(ctimeLsqnonlin)),ctimeSimP];
              [zeros(1,length(ctimeLsqnonlin)),ctimeSimQ]];
        description={'\\UseVerb{fmincon}';...
                     'Evo Sys $\\mathbb{P}$';...
                     'Evo Sys $\\mathbb{Q}$';...
                     'Simulation $\\mathbb{P}$';...
                     'Simulation $\\mathbb{Q}$'};
        file = fopen([tempPath,'\',latexFilePath],'a'); 
        fprintf(file,...
                '\\newcommand{%s}{\n',['\ctimeFminTotal',identifier,'value']);
        fprintf(file,'%3.3g',sum(ctimeLsqnonlin,'all'));
        fprintf(file,...
                '}\n');
        fprintf(file,...
                '\\newcommand{%s}{\n',['\ctimeEvoPTotal',identifier,'value']);
        fprintf(file,'%3.3g',ctimeEvoSysP);
        fprintf(file,...
                '}\n');
        fprintf(file,...
                '\\newcommand{%s}{\n',['\ctimeEvoQTotal',identifier,'value']);
        fprintf(file,'%3.3g',ctimeEvoSysQ);
        fprintf(file,...
                '}\n');
        fprintf(file,...
                '\\newcommand{%s}{\n',['\ctimeSimPTotal',identifier,'value']);
        fprintf(file,'%3.3g',ctimeSimP);
        fprintf(file,...
                '}\n');
        fprintf(file,...
                '\\newcommand{%s}{\n',['\ctimeSimQTotal',identifier,'value']);
        fprintf(file,'%3.3g',ctimeSimQ);
        fprintf(file,...
                '}\n');
        fprintf(file,...
                '\\SaveVerb{fmincon}=fmincon=\n');
        fprintf(file,...
                '\\newcommand{%s}{\n',latexCommand{1});
        fprintf(file,...
                '\\begin{tabular}{|c|*{%d}{c}|c|}\n',length(ratingMatrixYears));
        fprintf(file,'\\hline\n');
        % Header
        fprintf(file,'\\diagbox{Time}{Interval} & ');
        time=cat(2,0,ratingMatrixYears);
        for i=1:1:length(time)-1
            fprintf(file,'[%1.3f,%1.3f] & ',time(i),time(i+1));
        end
        fprintf(file,'Total\\\\\n');
        fprintf(file,'\\hline\n');
        
        % Body
        fprintf(file,mat2TableBody(description,data,'g',3,3,''));
            
        fprintf(file,'\\\\\\hline\n');
        fprintf(file,...
                '\\end{tabular}\n');
        fprintf(file,...
                '}\n');
        fclose(file);
    end
%% Rating Matrices
    function [latexFilePath,latexCommand]=saveRatingMatrices(saveAt)
        if strcmp(saveAt, '')
            latexFilePath='ratingMatrices.tex';
        else
            latexFilePath=[saveAt,'\','ratingMatrices.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        latexCommandTemp={['\ratingMatricesMarket',identifier];...
                      ['\ratingMatricesAdjusted',identifier];
                      ['\ratingMatricesAnalyticP',identifier];
                      ['\ratingMatricesAnalyticQ',identifier];
                      ['\ratingMatricesEvoP',identifier];
                      ['\ratingMatricesEvoQ',identifier];
                      ['\ratingMatricesSimP',identifier];
                      ['\ratingMatricesSimQ',identifier]};
        latexCommand={};
        description=ratings;
        ratingMatrices=zeros(size(Padjusted,1),size(Padjusted,2),size(Padjusted,3),8);
        ratingMatrices(:,:,:,1)=Pmarket;
        ratingMatrices(:,:,:,2)=Padjusted;
        ratingMatrices(:,:,:,3)=UPcal;
        ratingMatrices(:,:,:,4)=UQcal;
        ratingMatrices(:,:,:,5)=UP(:,:,tk);
        ratingMatrices(:,:,:,6)=UQ(:,:,tk);
        ratingMatrices(:,:,:,7)=XPProbs(:,:,tk);
        ratingMatrices(:,:,:,8)=XQProbs(:,:,tk);
        file = fopen([tempPath,'\',latexFilePath],'a'); 
        for r=1:1:length(latexCommandTemp)
            for k=1:1:size(Padjusted,3)
                latexCommand{end+1}=...
                    sprintf('%s%s',latexCommandTemp{r},numberDict{k});
                fprintf(file,...
                        '\\newcommand{%s%s}{\n',latexCommandTemp{r},numberDict{k});
                fprintf(file,...
                        '\\begin{tabular}{|c|*{%d}{c}|}\n',length(ratings));
                fprintf(file,'\\hline\n');
                % Header
                fprintf(file,'\\diagbox{From}{To} & ');
                for i=1:1:length(ratings)
                    if i<length(ratings)
                        fprintf(file,'%s & ',ratings{i});
                    else
                        fprintf(file,'%s\\\\ ',ratings{i});
                    end
                end
                fprintf(file,'\\hline\n');

                % Body
                fprintf(file,...
                    mat2TableBody(description,...
                                  ratingMatrices(:,:,k,r).*100,...
                                  'f',3,3,'\\,\\%%'));
                fprintf(file,'\\\\\\hline\n');
                fprintf(file,...
                        '\\end{tabular}\n');
                fprintf(file,...
                        '}\n');
                fprintf(file,...
                        '\\newcommand{%s%sCaption}{\n',latexCommandTemp{r},numberDict{k});
                fprintf(file,...
                    '\\caption{%s at %1.3f}',...
                    replace(replace(latexCommandTemp{r},'\ratingMatrices',''),identifier,''),...
                    ratingMatrixYears(k));
                fprintf(file,...
                        '}\n');
            end
        end
        fclose(file);
    end
%% Default Probabilities
    function [latexFilePath,latexCommand]=saveDefaultProbabilities(saveAt)
        if strcmp(saveAt, '')
            latexFilePath='defaultProbabilities.tex';
        else
            latexFilePath=[saveAt,'\','defaultProbabilities.tex'];
            mkDir([tempPath,'\',saveAt]);
        end
        latexCommandTemp={['\defaultProbabilities',identifier]};
        latexCommand={};
        description=ratings;
        file = fopen([tempPath,'\',latexFilePath],'a'); 
        for r=1:1:length(latexCommandTemp)
            for k=1:1:size(Padjusted,3)
                latexCommand{end+1}=...
                    sprintf('%s%s',latexCommandTemp{r},numberDict{k});
                fprintf(file,...
                        '\\newcommand{%s%s}{\n',latexCommandTemp{r},numberDict{k});
                fprintf(file,...
                        '\\begin{tabular}{|c|*{%d}{c}|}\n',1);
                fprintf(file,'\\hline\n');
                % Header
                fprintf(file,'\\diagbox{From}{To} & %s \\\\\n',ratings{end});
                fprintf(file,'\\hline\n');

                % Body
                fprintf(file,...
                    mat2TableBody(description,...
                                  PD(:,k).*100,...
                                  'f',3,3,'\\,\\%%'));
                fprintf(file,'\\\\\\hline\n');
                fprintf(file,...
                        '\\end{tabular}\n');
                fprintf(file,...
                        '}\n');
                fprintf(file,...
                        '\\newcommand{%s%sCaption}{\n',latexCommandTemp{r},numberDict{k});
                fprintf(file,...
                    '\\caption{%s at %1.3f}',...
                    replace(replace(latexCommandTemp{r},'\defaultProbabilities','Default Probabilities'),identifier,''),...
                    ratingMatrixYears(k));
                fprintf(file,...
                        '}\n');
            end
        end
        fclose(file);
    end
%% Figures
    function latexFilePath=saveFigures(figures,saveAt,figName)
        if strcmp(figName,'')
            figName='fig';
        end
        if strcmp(saveAt, '')
            latexFilePath='figure.tex';
            latexFolderPath=tempPath;
        else
            latexFilePath=[saveAt,'\','figure.tex'];
            latexFolderPath=[tempPath,'\',saveAt];
            mkDir([tempPath,'\',saveAt]);
        end
        file = fopen([tempPath,'\',latexFilePath],'a');
        for i=1:1:length(figures)
            picName=sprintf('%s_%d',figName,i);
            picPath = [latexFolderPath,'\',picName,'.',picType];
            if strcmp(saveAt, '')
                relPicPath=picName;
            else
                relPicPath=[saveAt,'\',picName];
            end
            saveas(figures(i),picPath,saveParam);
            fprintf(file,...
                '\\begin{landscape}\n');
            fprintf(file,...
                '\\includegraphics[width=.95\\columnwidth]{%s}\n',...
                changeSlash(relPicPath));
            fprintf(file,...
                '\\end{landscape}\n');
        end
        fclose(file);
    end
end
%% Auxiliary functions
function str=changeSlash(str)
    for i=1:1:length(str)
        if strcmp(str(i),'\')
            str(i)='/';
        end
    end
end
function delFile(file)
    if exist(file)
        delete(file);
    end
end
function delDir(dir)
    if exist(dir)==7
        rmdir(dir,'s');
    end
end
function cleanDir(mdir,except)
    except{end+1}='.';
    except{end+1}='..';
    for d = dir(mdir).'
      if ~any(strcmp(d.name,except))
          if d.isdir
              rmdir([d.folder,'\',d.name],'s');
          else
              delete([d.folder,'\',d.name]);
          end
      end
    end
end
function mkDir(dir)
    if exist(dir)==0
        mkdir(dir);
    end
end
function latexStr=mat2TableBody(description,mat,format,precision1,precision2,unit)
%     percent='\\\\,\\\\%%%%';
    latexStr='';
    for i=1:1:size(mat,1)
        if ~isempty(description)
            latexStr=[latexStr,description{i},' & '];
        end
        for j=1:1:size(mat,2)
            if j < size(mat,2)
                temp=sprintf(sprintf('$%%*.*%s$%%s & ',format),precision1,precision2,mat(i,j),unit);
            else
%                 temp=sprintf('$%*.*e$%s ',precision1,precision2,mat(i,j),unit);
                temp=sprintf(sprintf('$%%*.*%s$%%s ',format),precision1,precision2,mat(i,j),unit);
            end
            latexStr=[latexStr,temp];
        end
        if i < size(mat,1)
            temp='\\\\\n';
        else
            temp='\n';
        end
        latexStr=[latexStr,temp];
    end
end