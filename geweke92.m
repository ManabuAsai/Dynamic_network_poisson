function pval = geweke92(vX)
  iN = length(vX);
  iN1 = floor(0.1*iN); iN2 = floor(0.5*iN);
  vX1 = vX(1:iN1); vX2 = vX((iN-iN2+1):end);
  iX1bar = mean(vX1); iX2bar = mean(vX2);
  var1 = TSvar(vX1,0.1*iN1); var2 = TSvar(vX2,0.1*iN2);
  iZ = (iX1bar-iX2bar)/sqrt(var1/iN1 + var2/iN2);
  %iS2 = (iN1*var1 + iN2*var2)/(iN1+iN2);
  %iZ = (iX1bar-iX2bar)/sqrt(iS2);
  pval = 2*(1-cdf('norm',abs(iZ),0,1));
end