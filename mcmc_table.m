function mO = mcmc_table(mR)
[iN,iP]=size(mR);
mRR = sort(mR);
mO = [mean(mRR)' mRR(ceil(0.025*iN),:)' mRR(floor(0.975*iN),:)']; 
vG = zeros(iP,2);
for i=1:iP
    vG(i,1) = geweke92(mR(:,i));
    vG(i,2) = TSvar(mR(:,i));
end
mO = [mO vG(:,1) (var(mR)'./vG(:,2))];
end