function [vTh_new,iA] = DRAM_dnar1wi(vTheta,mTheta,j,mS1,mY,mDY,mWDY,iPar,vMu0,mLam0,vLam_ini)
    %iNs = 0.5*iN*(iN+1);
    iP = length(vTheta);sp=(2.4)^2/iP;
    eps = min(diag(cov(mTheta)));
    vD = vTheta;
    if (j<1500)
        mS = sp*mS1;
    else
        mS = sp*(cov(mTheta)+0.01*eps*eye(iP));
    end
    mS1h = real(sqrtm(mS));
    iY0 = -0.5*log(det(mS))-0.5*(vTheta-vD)'*(mS\(vTheta-vD));
    gam0=0.01;
    iY02 = -0.5*log(det(gam0*mS))-0.5*(vTheta-vD)'*((gam0*mS)\(vTheta-vD));
    iX0 = loglike_dnar1wi(mY,mDY,mWDY,vTheta,iPar,vLam_ini)+LogPrior(vTheta,vMu0,mLam0);
    iG0 = iX0-iY0;
    iI=100;
    vA = zeros(iI,1); vGc = zeros(iI,1);
    for i=1:iI
        if (i==1)
            gam=1;
        else
            gam = gam0; iG0 = iX0-iY02;
        end
        vTh_c = vD + sqrt(gam)*mS1h*random('Normal',0,1,iP,1);
        ll = loglike_dnar1wi(mY,mDY,mWDY,vTh_c,iPar,vLam_ini);
        %h=0;
        %while (h<1)
        %    vTh_c = vD + sqrt(gam)*mS1h*random('Normal',0,1,iP,1);
        %    %vTh_c(iP)
        %    iZ = prod(abs(vTh_c(iNs+1:2*iNs))<1);
        %    ll = loglike_dnar1w(mY,mDY,mWDY,vTh_c,iPar);
        %    if (iZ==1)&&(ll>-1e+10)
        %        h=1;
        %    end
        %end
        %vTh_c
        iYc = -0.5*log(det(gam*mS))-0.5*(vTh_c-vD)'*((gam*mS)\(vTh_c-vD));
        iXc = LogPrior(vTh_c,vMu0,mLam0) +ll;
        iGc = iXc-iYc;
        vGc(i) = exp(iGc-iG0);%[iGc-iG0]
        if (i==1)
            alpha_mh = min([1 vGc(i)]); vA(1) = alpha_mh;
        else
            alpha_mh = min([1 max([0 (vGc(i) - max(vGc(1:i-1)))])/(1-max(vGc(1:i-1)))]);
            vA(i) = alpha_mh;
        end
        %[vGc(i) alpha_mh]
        u = rand; iA = (u<alpha_mh);%[i alpha_mh]
        if (iA==1)
            break;
        end
    end
    if (iA==1)
        vTh_new = vTh_c;
    else
        vTh_new = vTheta;
    end
end
    
        
