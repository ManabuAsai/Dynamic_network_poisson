function vH = vech(mX)
  n = length(mX);
  b =0;
  vH = NaN(n*(n+1)/2,1);
  vH(1:n) = mX(1:n);
  idx=n+1;
  for j=2:n;
      vH(idx:idx+(n-j))=mX(j:end,j);
      idx = idx + n-j+1;
  end
end