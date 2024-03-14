function mY = DGP_TNW(mDN,gam,alpha,beta,psi,i)
% DGP based on the true network weight
% p=1
  rng(2*i+123,'twister');
  [iT,iN2]=size(mDN); iN=sqrt(iN2);
  mY=zeros(iT,iN);
  vYst = zeros(1,iN); vWYst = zeros(1,iN); vLam = zeros(1,iN);
  %vYst = random('Poisson',10*ones(1,iN));
  for t=1:iT
      vLam = alpha*vLam + gam*ones(1,iN) + beta*vYst + psi*vWYst;
      vY = random('Poisson',exp(vLam));
      mY(t,:) = vY; vYst = log(1+vY);
      mGt = reshape(mDN(t,:)',iN,iN);
      vGt = sum(mGt,2); vGt = (vGt < 1) + vGt;
      mWt = inv(diag(vGt))*mGt;
      vWYst = vYst*mWt'; 
  end
end
     
  