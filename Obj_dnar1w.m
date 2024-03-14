function obj = Obj_dnar1w(vP,mY,mDY,mWDY,iP,vLam_ini)
  [iT,iN] = size(mY);
  ll = loglike_dnar1w(mY,mDY,mWDY,vP,iP,vLam_ini);
  obj = - ll/((iT*iN));%obj
end