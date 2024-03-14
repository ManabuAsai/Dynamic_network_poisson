% Calculation of Top Lyapunov exponent

% upload linkage variables
load('linkvar1.mat');% mDN0
mDN = mDN0(1:end-30,:);

[iT,iN0]=size(mDN);iN=51;
mLy_exp = zeros(iT,3);
for i=1:3
    iI=i;
    if (iI==1)
       gam = 0.6; beta = 0.93; psi = -0.01; alpha = 0.02;
    elseif (iI==2)
       gam = 0.6; beta = 0.90; psi = -0.01; alpha = 0.07;
    elseif (iI==3)
       gam = 0.6; beta = 0.93; psi = -0.03; alpha = 0.02;
    end
    vLy_exp = zeros(iT,1);
    mA = alpha*eye(iN);mAm=eye(iN);
    for t=1:iT
        mGt = reshape(mDN(t,:)',iN,iN); mGt = unvech(vech(mGt))+unvech(vech(mGt'));
        vGt = sum(mGt,2); vGt = (vGt < 1) + vGt;
        mWt = inv(diag(vGt))*mGt;
        mCt = beta*eye(iN)+psi*mWt;
        mAm = (abs(mA)+abs(mCt))*mAm;
        iRho = sqrt(max(eig(mAm'*mAm)));
        vLy_exp(t) = (1/t)*log(iRho);
    end
    mLy_exp(:,i)=vLy_exp;
end
plot((1:iT-1),mLy_exp(2:end,1:3));xlim([0,80]);ylim([-0.06,0]);
legend("DGP1","DGP2","DGP3",'Location','southeast')


