% data transformation
%
clear
clc
iK=2;
if (iK==1)
    load('linkvar.mat');% mD
    mDN0 = mD(1:end-6,:)+mD(2:end-5,:)+mD(3:end-4,:)+mD(4:end-3,:)+mD(5:end-2,:)+mD(6:end-1,:)+mD(7:end,:);
    mDN0 = mDN0(4:7:end,:); % sum at Sunday to Saturday
    mDN = mDN0(end-81-30+1:end,:);
    mDN0 = (mDN>0);
    save('linkvar1.mat','mDN0');
elseif (iK==2)
    load('confirmed.mat');% mD1
    mY0 = mD1(1:end-6,:)+mD1(2:end-5,:)+mD1(3:end-4,:)+mD1(4:end-3,:)+mD1(5:end-2,:)+mD1(6:end-1,:)+mD1(7:end,:);
    mY0 = mY0(4:7:end,:); % sum at Sunday to Saturday
    mY0 = mY0(end-81-30+1:end,:);
    save('confirmed1.mat','mY0');
end

