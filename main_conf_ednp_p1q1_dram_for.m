% EN-PAR(1,1)
% weekly data
% Estimation via DRAM
% One-step-ahead forecast
clear
clc
load('confirmed1.mat');
mY = mY0;
mLY = log(1+mY);

[iT0,iN] = size(mY);
iP=1;
iT = iT0-iP;
[iT iN]

iM=6;
mDN = DyNet(mLY,iM);
mWLY = zeros(iT0,iN);
for t=1:iT0
   mGt = unvech(mDN(t,:)');
   mGt = (mGt<0)+mGt;
   vGt = sum(mGt,2); vGt = (vGt < 1) + vGt;
   mWt = inv(diag(vGt))*mGt;
   vWYt = mLY(t,:)*mWt'; 
   mWLY(t,:) = vWYt;
end

vP0= zeros(2*iP+2,1); mCP0 = 5*eye(2*iP+2);
mC0=0.25*eye(2*iP+1);
rng(123,'twister');

%options = optimset('Display','none','Algorithm','interior-point','MaxIter',1e5,'MaxFunEvals',1e5);
options = optimset('Display','iter','Algorithm','interior-point','MaxIter',1e5,'MaxFunEvals',1e5);
lb = -5*ones(2*iP+2,1);
ub =  5*ones(2*iP+2,1);
vP_ini = [0.1; 1e-06*ones(2*iP,1);0.1];

st=500;iL=1000;
iR= st+iL;%iR=10;st=0;

iT=81-iP;iF=30;
mFFm = zeros(iF,iN);mFFl = zeros(iF,iN);mFFu = zeros(iF,iN);
for k=1:iF
    mY0 = mY(end-iT-iF+k-iP:end-iF+k-1,:);
    mLY0 = mLY(end-iT-iF+k-iP:end-iF+k-1,:);
    
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
    
    mWLY0 = mWLY(end-iT-iF+k-iP:end-iF+k-1,:);
    vLam_ini = 1 + mY(end-iT-iF+k-1,:);

    [vP_h,~] = fmincon(@Obj_dnar1w,vP_ini,[],[],[],[],lb,ub,[],options,mY0,mLY0,mWLY0,iP,vLam_ini);
    mCh = -inv(D2_dnar1w(vP_h,mY0,mLY0,mWLY0,iP,vLam_ini));
    [mP,mE] = eig(mCh);
    if (min(diag(mE))<0)
        mCh = mP*abs(mE)*mP';
    end

    %vP = vP_ini;
    vP = vP_h;% for smaller iL 
    mS1=mCh;
    mF = zeros(iR,iN);
    mRR = zeros(iR,2*iP+2);
    alpha_mh =0;
    for i=1:iR
        [vP,a_mh] = DRAM_dnar1w(vP,mRR,i,mS1,mY0,mLY0,mWLY0,iP,vP0,mCP0,vLam_ini);
        mRR(i,:)=vP';
        alpha_mh = alpha_mh + a_mh;
        vLLam = lambda_dnar1w_1for(mY0,mLY0,mWLY0,vP,iP,vLam_ini);
        mF(i,:) = exp(vLLam)-1;
    end
    alpha_mh = alpha_mh/iR;alpha_mh
    mY_for = sort(mF(st+1:end,:));
    mFFm(k,:) = mY_for(0.5*iL,:);
    mFFl(k,:) = mY_for(0.05*iL,:);
    mFFu(k,:) = mY_for(0.95*iL,:);
end
save('mR_conf_ednp_p1q1_for','mFFm',"mFFu","mFFl");

