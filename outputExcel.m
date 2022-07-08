function outputExcel(fileName,...
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
                   preDefaultPlots,...
                   r,m,LGDI,LGDC,...
                   thresholdsIUC,thresholdsCUC,thresholdsIRT,thresholdsCRT,thresholdsIPC,thresholdsCPC,...
                   collateralPlotsUC,collateralPlotsRT,collateralPlotsPC,...
                   cvaPlotsUC,cvaPlotsRT,cvaPlotsPC)
   

root=[pwd, '\' ,'Results'];
excelRoot=[root,'\','Excel'];
dataPath=[excelRoot,'\',fileName];
fileName=[dataPath,'\',fileName,'.xlsx'];
          
delDir(dataPath)
mkDir(dataPath)

matrix3dExcel('Generator P',AP);
matrix3dExcel('Generator Q',AQ);
matrix3dExcel('Rating Matrices Analytic P',UPcal);
matrix3dExcel('Rating Matrices Analytic Q',UQcal);
matrix3dExcel('Rating Matrices Simulated P',XPProbs(:,:,tk));
matrix3dExcel('Rating Matrices Simulated Q',XQProbs(:,:,tk));
matrix3dExcel('Rating Matrices Euler P',UP(:,:,tk));
matrix3dExcel('Rating Matrices Euler Q',UQ(:,:,tk));
matrix3dExcel('Adjusted Rating Matrices P',Padjusted);
matrix3dExcel('Market Rating Matrices P',Pmarket);
default2Excel('Market Default Probabilities Q',PD);

    function matrix3dExcel(sheet,data)
        for i=1:1:size(data,3)
            description={sprintf('At time %3.3g',ratingMatrixYears(i))};
            T=array2table(data(:,:,i),'VariableNames',ratings,...
                                      'RowNames',ratings,...
                                      'DimensionNames',{'From/To','Variables'});
            writetable(T,fileName,'Sheet',sheet,...
                                  'WriteMode','append',...
                                  'WriteRowNames',true,...
                                  'WriteVariableNames',true);
            writecell(description,fileName,'Sheet',sheet,'WriteMode','append');
        end
    end
    function default2Excel(sheet,data)
        for i=1:1:size(data,2)
            description={sprintf('At time %3.3g',ratingMatrixYears(i))};
            T=array2table(data(:,i),'VariableNames',ratings(end),...
                                      'RowNames',ratings,...
                                      'DimensionNames',{'From/To','Variables'});
            writetable(T,fileName,'Sheet',sheet,...
                                  'WriteMode','append',...
                                  'WriteRowNames',true,...
                                  'WriteVariableNames',true);
            writecell(description,fileName,'Sheet',sheet,'WriteMode','append');
        end
    end
   
end

%% Auxiliary functions
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