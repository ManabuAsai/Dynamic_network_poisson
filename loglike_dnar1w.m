function ll = loglike_dnar1w(mY,mLY,mWLY,vP,iP,vLam_ini)
  % no indivual effects
  [iT0,iN]=size(mY); iT=iT0-iP;
  if nargin < 6, vLam_ini = ones(1,iN); end
  vB = vP(1:end-1);
  phi = vP(end);

  loglike = 0;mLLam=zeros(iT+1,iN);mLLam(1,:) = log(vLam_ini);
  for t=1:iT
    mYYt = zeros(iN,iP); mWWt = zeros(iN,iP);
    for i=1:iP
        mYYt(:,i) = mLY(iP+t-i,:)';
        mWWt(:,i) = mWLY(iP+t-i,:)';
    end
    mXt = [ones(iN,1) mYYt mWWt];
    logLambda = phi*mLLam(t,:) + vB'*mXt';
    mLLam(1+t,:) = logLambda;
    loglike = loglike + sum(mY(iP+t,:).*logLambda - exp(logLambda)-gammaln(mY(iP+t,:)+1));
    %loglike = loglike + sum(mY(iP+t,:).*logLambda - exp(logLambda));
  end
  ll = loglike;
end
