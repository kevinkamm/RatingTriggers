function r=fminJLT(h,Amarket,PD,t)
h=[reshape(h,[],1);1];
Ah=diag(h)*Amarket;
leftQ=expm(Ah.*t);
rightQ=PD;

r=sum(abs(leftQ(1:end-1,end)-rightQ(1:end-1)).^2);
end