% CN-PAR model
% DRAM estimation
% Monte Carlo experiment
% Use constant network when the DGP is generated via the true dynamic network

clear
clc
iT=80; iN=51;
p=0.6;

iT0=iT+1;

iP=1;
%gam = 0.6;beta = 0.93;psi = -0.01;alpha = 0.02;% DGP1
%gam = 0.6;beta = 0.9;psi = -0.01;alpha = 0.07;% DGP2
gam = 0.6;beta = 0.93;psi = -0.03;alpha = 0.02;% DGP3

vTheta0=[gam;beta;psi;alpha];

load('linkvar1.mat');
mTDN = mDN0(1:end-30,:);
for t=1:iT0
   mGt = reshape(mDN1(t,:)',iN,iN); mGt = unvech(vech(mGt))+unvech(vech(mGt'));
   mTDN(t,:) = vec(mGt)';
end

iQ=500;
mQ = zeros(iQ,2*iP+2);
for q=1:iQ
    q
mY = DGP_TNW(mTDN,gam,alpha,beta,psi,q);
mLY = log(1+mY);

mDY = mY(2:end,:)-mY(1:end-1,:);
mCorr = corrcov(cov(mDY));
iC = tinv(0.95,length(mDY)-2);
mC2=zeros(iN,iN);
for i=2:iN
    for j=1:(i-1)
        iRho=mCorr(i,j);
        mC2(i,j)=(abs(iRho*sqrt((length(mDY)-2)/(1-iRho^2)))>iC);
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
lb = -5*ones(2*iP+2,1);
ub =  5*ones(2*iP+2,1);
vP_ini = [gam;beta;psi;alpha];
vLam_ini = 1+mY(iP,:);
[vP_h,fval] = fmincon(@Obj_dnar1w,vP_ini,[],[],[],[],lb,ub,[],options,mY,mLY,mWLY,iP,vLam_ini);
format longG
mNH = -(1/iT)*D2_dnar1w(vP_h,mY,mLY,mWLY,iP,vLam_ini);
if (min(eig(mNH))<0)
    mNH = 0.8*mNH + 0.2*diag(diag(mNH));
end
mCh = inv(mNH);
[vP_h sqrt(diag(mCh))]

vP0= zeros(2*iP+2,1); mCP0 = 5*eye(2*iP+2);
rng(123,'twister');

st=10000;
iR= st+20000;
alpha_mh=0;j=0;
mRR = zeros(iR,1+2*iP+1);
mS1=1*mCh; vP=vP_h;
for i=1:iR
  j=j+1;
  [vP,a_mh] = DRAM_dnar1w(vP,mRR,j,mS1,mY,mLY,mWLY,iP,vP0,mCP0,vLam_ini);
  mRR(i,:) = vP';
  alpha_mh = alpha_mh + a_mh;
  if (j==1000)
    i
    if (i<st+10)
        [mean(mRR(1:i,2:iP+1))'     mean(mRR(1:i,iP+2:2*iP+1))']
        [mean(mRR(1:i,1))     mean(mRR(1:i,end))]
    else
        [mean(mRR(st+1:i,2:iP+1))'     mean(mRR(st+1:i,iP+2:2*iP+1))']
        [mean(mRR(st+1:i,1))     mean(mRR(st+1:i,end))]
    end
    j=0;
  end

end
alpha_mh=alpha_mh/iR;alpha_mh
mR = mRR(st+1:end,:);
mQ(q,:)=mean(mR);
end
[mean(mQ); sqrt(var(mQ)); sqrt(mean((mQ-ones(iQ,1)*vTheta0').^2))]'

