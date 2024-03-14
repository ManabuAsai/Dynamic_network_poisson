function vLLam = lambda_dnar1wi_1for(mY,mDY,mWDY,vP,iP,vLam_ini)
  % no indivual effects
  [iT0,iN]=size(mY); iT=iT0-iP;
  if nargin < 6, vLam_ini = ones(7,iN); end
  vB = vP(1:end-1);
  phi = vP(end);

  mLLam=zeros(iT+2,iN);mLLam(1,:) = log(vLam_ini);
  for t=1:iT+1
    mYYt = zeros(iN,iP); mWWt = zeros(iN,iP);
    for i=1:iP
        mYYt(:,i) = mDY(iP+t-i,:)';
        mWWt(:,i) = mWDY(iP+t-i,:)';
    end
    mXt = [eye(iN) mYYt mWWt];
    logLambda = phi*mLLam(t,:) + vB'*mXt';
    mLLam(1+t,:) = logLambda;
  end
vLLam = mLLam(end,:);
end
