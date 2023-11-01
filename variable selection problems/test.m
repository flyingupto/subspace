%%%%%%%%%%%%08-06-2023 by Wang Shangfei%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all;
load('dataset.mat');
nsamples=10000;
burnin=5000;
thin=5;
[x1,x2,time] = bhs_compare(A, y_noisy,nsamples,burnin, thin);
meanx1=mean(x1,2);
meanx2=mean(x2,2);
[x,fitinfo] = lasso(A,y_noisy,'CV',10); 
lam = fitinfo.IndexMinMSE;  
mat = x(:,lam);             
[row, ] = find(x(:,lam)~=0);

[xmax1,I1]=sort(abs(meanx1));
num=9;
sn=sqrt(num);
xmax1=xmax1(end-num+1:end);
I1=I1(end-num+1:end);

data=[meanx1 meanx2 mat];
figure;
for i=1:num
subplot(sn,sn*2,2*i-1);
[f1,xi1] = ksdensity(x1(I1(i),:));
plot(xi1,f1,'Color', [0 0.4470 0.7410]);
title(I1(i));
hold on
fm=max(f1);
plot([meanx2(I1(i)),meanx2(I1(i))],[0,fm],'k');
hold on
plot([meanx1(I1(i)),meanx1(I1(i))],[0,fm],'Color', [0 0.4470 0.7410]);
hold on
plot([mat(I1(i)),mat(I1(i))],[0,fm],'r');
hold on


subplot(sn,sn*2,2*i);
[f2,xi2] = ksdensity(x2(I1(i),:));
plot(xi2,f2,'k');
title(I1(i));
hold on
fm=max(f2);
plot([meanx1(I1(i)),meanx1(I1(i))],[0,fm],'Color', [0 0.4470 0.7410]);
hold on
plot([mat(I1(i)),mat(I1(i))],[0,fm],'r');
hold on
plot([meanx2(I1(i)),meanx2(I1(i))],[0,fm],'k');
hold on
end

figure;
for i=1:num
subplot(sn,sn*2,2*i-1);
plot(1:size(x1,2),x1(I1(i),:),'Color', [0 0.4470 0.7410]);
title(I1(i));
hold on;
subplot(sn,sn*2,2*i);
plot(1:size(x2,2),x2(I1(i),:),'k');
title(I1(i));
hold on;
end