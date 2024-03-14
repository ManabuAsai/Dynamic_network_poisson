function  mD2 = D2_dnar1w(psi,mY,mDY,mWDY,iPar,vLam_ini)
  iP = length(psi);
  h = abs(psi)*10^(-4);
  mD2 = zeros(iP,iP);
  for i=1:iP
      for j=1:i
      psi_ij = psi;psi_i = psi; psi_j = psi;  
      psi_ij(i) = psi(i) + h(i); psi_ij(j) = psi_ij(j)+h(j);
      psi_i(i) = psi(i)+h(i);
      psi_j(j) = psi(j)+h(j);
      c1 = loglike_dnar1w(mY,mDY,mWDY,psi_ij,iPar,vLam_ini);
      c2 = loglike_dnar1w(mY,mDY,mWDY,psi_i,iPar,vLam_ini);
      c3 = loglike_dnar1w(mY,mDY,mWDY,psi_j,iPar,vLam_ini);
      c4 = loglike_dnar1w(mY,mDY,mWDY,psi,iPar,vLam_ini);
      mD2(i,j) = (1/(h(i)*h(j)))*(c1-c2-c3+c4);
      end
  end
  mD2 = unvech(vech(mD2));
end
       