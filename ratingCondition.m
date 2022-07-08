function c=ratingCondition(A)
K=size(A,1);
c=zeros(K-1,K-1);
for i=1:1:K-1
    for k=1:1:K-1
        if k==i+1
            continue;
        end
        s1=sum(A(i,k:end),2);
        s2=sum(A(i+1,k:end),2);
        c(i,k)=s1-s2;
    end
end
% diagInd=eye(size(A,1),'logical');
% A(diagInd)=0;
% Ac=flip(cumsum(flip(A,2),2),2);
% c=Ac(1:end-1,:)-Ac(2:end,:);
% K=size(A,1);
% for i=1:1:K-1
%     for k=1:1:K
%         s1=0;
%         for j=k:1:K
%             if j==i
%                 continue;
%             end
%             s1=s1+A(i,j);
%         end
%         s2=0;
%         for j=k:1:K
%             if j==i+1
%                 continue;
%             end
%             s2=s2+A(i+1,j);
%         end
%         c(i,k)=s1-s2;
%     end
% end
end