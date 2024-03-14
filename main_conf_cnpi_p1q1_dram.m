% CN-PAR-I(1,1)
% weekly data
% Estimation via DRAM
% individual effects
clear
clc
load('confirmed1.mat');
mY = mY0(1:end-30,:);
mLY = log(1+mY);
iP=1;

[iT0,iN] = size(mY);[iT0,iN]

mCorr = corrcov(cov(mLY));
iC = tinv(0.95,length(mLY)-2);
mC2=zeros(iN,iN);
for i=2:iN
    for j=1:(i-1)
        iRho=mCorr(i,j);
        mC2(i,j)=(abs(iRho*sqrt((length(mLY)-2)/(1-iRho^2)))>iC);
    end
end
mC2 = unvech(vech(mC2));
mDN = ones(iT0,1)*vec(mC2)';

mWLY = zeros(iT0,iN);
for t=1:iT0
   mGt = reshape(mDN(t,:)',iN,iN);
   vGt = sum(mGt,2); vGt = (vGt < 1) + vGt;
   mWt = inv(diag(vGt))*mGt;
   vWYt = mLY(t,:)*mWt'; 
   mWLY(t,:) = vWYt;
end


options = optimset('Display','iter','Algorithm','interior-point','MaxIter',1e5,'MaxFunEvals',1e5);
lb = -5*ones(2*iP+1+iN,1);
ub =  5*ones(2*iP+1+iN,1);
vP_ini = [0.5*ones(iN,1); 1e-05*ones(2*iP,1);-0.1];
vLam_ini = 1+mY(iP,:);
[vP_h,fval] = fmincon(@Obj_dnar1wi,vP_ini,[],[],[],[],lb,ub,[],options,mY,mLY,mWLY,iP,vLam_ini);
format longG
mCh = -inv(D2_dnar1wi(vP_h,mY,mLY,mWLY,iP,vLam_ini));
[vP_h sqrt(diag(mCh))]

vP0= zeros(2*iP+1+iN,1); mCP0 = 5*eye(2*iP+1+iN);
rng(123,'twister');


st=20000;
iR= st+20000;
alpha_mh=0;j=0;
mRR = zeros(iR,iN+2*iP+1);
mS1=mCh; vP=vP_h;
for i=1:iR
  j=j+1;
  [vP,a_mh] = DRAM_dnar1wi(vP,mRR,j,mS1,mY,mLY,mWLY,iP,vP0,mCP0,vLam_ini);
  mRR(i,:) = vP';
  alpha_mh = alpha_mh + a_mh;
  if (j==1000)
    i
    if (i<st+10)
        [mean(mRR(1:i,iN+1:iN+iP))'     mean(mRR(1:i,iN+iP+1:iN+2*iP))']
        [mean(mRR(1:i,1))     mean(mRR(1:i,end))]
    else
        [mean(mRR(st+1:i,iN+1:iN+iP))'     mean(mRR(st+1:i,iN+iP+1:iN+2*iP))']
        [mean(mRR(st+1:i,1))     mean(mRR(st+1:i,end))]
    end
    j=0;
  end

end
alpha_mh=alpha_mh/iR;alpha_mh
mR = mRR(st+1:end,:);
format long
mO = mcmc_table(mR);mO
dic = DIC_pr_dnar1wi(mY,mLY,mWLY,mR,iP,vLam_ini);dic

