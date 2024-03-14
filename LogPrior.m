function log_prior = LogPrior(vB,vB0,mCB0)
  iP = length(vB);
  log_prior = -0.5*iP*log(2*pi) -0.5*log(det(mCB0)) - 0.5*(vB-vB0)'*(mCB0\(vB-vB0));
end