function [Padjusted,A,Pmarket,ratings]=ratingMatrixLoader(ratingMatrixFolder,...
                                    ratingAgency,...
                                    ratingYear,...
                                    ratingMatrixDataset,...
                                    ratingMatrixYears)
%%RATINGMATRIXLOADER loads rating matrix from destination folder with 
% given data set and array of years
%   Input:
%       ratingMatrixFolder (str): contains the path to the folder of the rating
%                          matrices
%       ratingAgencey (str): contains the name of the rating agency
%       ratingYear (str): contains the year of the dataset
%       ratingMatrixDataset (int): contains the number of the data set
%       ratingMatrixYears (1xp array): contains the years of the rating matrices
%   Output:
%       Padjusted (KxKxp array): contains the rating probabilities in decimals
%       A (KxKxp array): contains the generator in decimals
%       Pmarket (KxKxp array): contains the rating probabilities in
%                              decimals without withdrawal adjustment
%       ratings (1xK cell array): contains the names of the ratings
Padjusted=[];
Pmarket=[];
A=[];
for i=1:1:length(ratingMatrixYears)
    ratingMatrixName=sprintf('%s_%s_%d_%2.2f*y.*',...
                             ratingAgency,...
                             ratingYear,...
                             ratingMatrixDataset,...
                             ratingMatrixYears(i));
    d=dir([pwd,'/',ratingMatrixFolder,'/',ratingAgency,'/',ratingMatrixName]);
    table=readtable([ratingMatrixFolder,'/',ratingAgency,'/',d.name],...
                             'VariableNamingRule','preserve');
    headers = table.Properties.VariableNames;
    data = table2array(table(1:end,2:end));
    if strcmp(headers{1},'%')
        data=data./100;
    end
    Pmarket(:,:,i)=data;
    Padjusted(:,:,i)=adjustmentWithdrawal(data);
    temp=logm(data);
    diagInd=eye(size(temp,1),'logical');
    temp(temp<=0)=0;
    temp(diagInd)=0;
    temp(diagInd)=-sum(temp,2);
    A(:,:,i)=temp;
end
ratings=headers(2:end);

end

%% Choose a rating matrix adjustment by (un-)commenting

% function MA=adjustmentWithdrawal(M)
% % adjustment matrix
% weightMatrix=full(spdiags([12.5,62.5,25].*ones(size(M,2),1),[-1 0 1],size(M,2),size(M,2)));
% weightMatrix(end,:)=0;
% weightMatrix(:,end)=0;
% weightMatrix(1,3)=12.5;
% weightMatrix(end-1,end-2)=25;
% weightMatrix(end-1,end-3)=12.5;
% weightMatrix=weightMatrix./100;
% wd=sum(M,2);
% MA=M+weightMatrix.*(1-wd);
% end


% function MA=adjustmentWithdrawal(M)
% % add not rated to diagonal only
% MA=M;
% wd=sum(M,2);
% for i=1:1:length(wd)
%     if wd(i)>0
%         MA(i,i)=M(i,i)+1-wd(i);
%     end
% end
% end

% function MA=adjustmentWithdrawal(M)
% % exclude default for adjustment
% MA=zeros(size(M));
% wd=sum(M,2);
% for i=1:1:length(wd)
%     if wd(i)>0
%         y=M(i,1:end-1);
%         y(y==0)=1e-10;
%         b=(y/sum(y,'all')).*(1-wd(i));
%         MA(i,1:end-1)=M(i,1:end-1)+b;
%     else
%         MA(i,:)=M(i,:);
%     end
% end
% MA(:,end)=M(:,end);
% end



function MA=adjustmentWithdrawal(M)
%the usual adjustment
MA=zeros(size(M));
wd=sum(M,2);
for i=1:1:length(wd)
    if wd(i)>0
        y=M(i,:);
        y(y==0)=1e-10;
        b=(y/sum(y,'all')).*(1-wd(i));
        MA(i,:)=M(i,:)+b;
    else
        MA(i,:)=M(i,:);
    end
end
end




% function MA=adjustmentWithdrawal(M)
% MA=zeros(size(M));
% wd=sum(M,2);
% for i=1:1:length(wd)
%     if wd(i)>0
%         x=1:1:length(wd);
%         y=M(i,:);
%         y(y==0)=1e-10;
%         deg=4;
%         lreg=polyfit(x,log(y),deg);
% %         temp=zeros(size(wd))';
% %         for j=1:1:deg
% %             temp=temp+lreg(j).*x.^(deg+1-j);
% %         end
% %         yreg=exp(lreg(end)).*exp(temp);
%         yreg=exp(polyval(lreg,x));
% %         figure();hold on;
% %         xx=linspace(1,length(wd),1000);
% %         yy=polyval([lreg(1:end-1),0],xx);
% %         plot(xx,exp(lreg(end)).*exp(yy)./...
% %             (((length(wd)-1)/(1001)).*sum(exp(lreg(end)).*exp(yy),'all')))
% %         plot(x,yreg./sum(yreg,'all'))
% %         plot(x,M(i,:)./sum(M(i,:),'all'))
%         b=(yreg/sum(yreg,'all')).*(1-wd(i));
%         MA(i,:)=M(i,:)+b;
%     else
%         MA(i,:)=M(i,:);
%     end
% end