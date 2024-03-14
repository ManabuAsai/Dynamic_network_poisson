function dic = DIC_pr_dnar1wi(mY,mDY,mWDY,mR,iP,vLam_ini)
  vP = mean(mR)';
  iD1 = loglike_dnar1wi(mY,mDY,mWDY,vP,iP,vLam_ini);
  iR = length(mR);
  vD2 = zeros(iR,1);
  for i=1:iR
      vD2(i) = loglike_dnar1wi(mY,mDY,mWDY,mR(i,:)',iP,vLam_ini);
  end
  vPd= -2*vD2+2*iD1;
  vDIC = 2*iD1 + vPd;
  [mean(vDIC) sqrt(var(vDIC)) max(vDIC) min(vDIC)]
  dic = mean(vDIC);
end
