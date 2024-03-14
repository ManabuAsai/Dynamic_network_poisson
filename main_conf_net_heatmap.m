% Heatmap for extracted network structure

str = {'AL','AK','AZ','AR','CA','CO','CT','DE','DC','FL','GA','HI','ID','IL','IN','IA','KS','KY','LA','ME','MD','MA','MI','MN','MS','MO','MT','NE','NV','NH','NJ','NM','NY','NC','ND','OH','OK','OR','PA','RI','SC','SD','TN','TX','UT','VT','VA','WA','WV','WI','WY'};


load('linkvar1.mat');
mDN = mDN0;

[iT0,~] = size(mDN);iN=51;
mDNT = zeros(iT0,iN*iN);
for t=1:iT0
   mDNt = reshape(mDN(t,:)',iN,iN); 
   mDNt = unvech(vech(mDNt))+unvech(vech(mDNt'));
   mDNT(t,:) = vec(mDNt)';
end

mADN =reshape(mean(mDNT),iN,iN);
h = heatmap(str,str,mADN);
xlabel('From:')
ylabel('To:')

