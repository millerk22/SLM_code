% normalized mutual information, as used in Lancichinette's benchmark
% paper.
% Take input C and Cp, both are N by 1 vector indicating the group
% assignment. 
% The first two lines are placed here To make sure C and C' take values 
% among {1,2,..} starting from 1 without skipping any numbers.

function NMI = nmiMeas(C,Cp)
	[~,~,C] = unique(C);
	[~,~,Cp] = unique(Cp);

N=size(C,1);
Nc=max(C);
Ncp=max(Cp);
Nprod=Nc*Ncp;

if Nc == 1 && Ncp == 1
	NMI = 1;
	return;
end

Pkp=hist(Cp,1:max(Cp))/N;

num = 0;
for k=1:max(C)
  index=find(C==k);
  Pk=size(index,1)/N;
  Pkkp=hist(Cp(index),1:max(Cp));
  for kp=find(Pkkp)
    Pkkp(kp)=Pkkp(kp)/N;
    num = num + Pkkp(kp)*log(Pkkp(kp)/(Pk*Pkp(kp)));
  end
end
num = 2*num;

denom = 0;
Pk =hist(C ,1:max(C ))/N;
for k=1:max(C)
    denom = denom + Pk(k)*log(Pk(k));
end
for kp=1:max(Cp)
    denom = denom + Pkp(kp)*log(Pkp(kp));
end
denom = -1*denom;

NMI = num/denom;
