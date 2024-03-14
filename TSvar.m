function var = TSvar(vX0,bm)
% compute spectral density function using bandwidth=bm at 2 points (0,pi)
% and Parzen window
iN = length(vX0); 
vX = vX0 - mean(vX0);
%sp = (1/iN)*periodogram(vX,bm); % Need Signal Processing Toolbox 
%var = 2*pi*sp(1);
  xdft = fft(vX);
  iL = floor(iN/2+1);
  xdft = xdft(1:iL,:);
  psdx = (1/(2*pi*iN))*(abs(xdft).^2);
  psdx(2:iL-1) = 2*psdx(2:iL-1);
  var = (1/iL)*sum((ones(iL,1)-(1/iL)*(0:iL-1)').*psdx);
end