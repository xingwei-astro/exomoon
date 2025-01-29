clear; close all; clc;

for i=1:100
map(i,1)=0.7; map(i,2)=0.7; map(i,3)=0.7;
end

load rTrH.dat
ni=101; nj=101;
T=reshape(rTrH(:,1),nj,ni);
Bp=reshape(rTrH(:,2),ni,nj);
rH=reshape(rTrH(:,3),ni,nj);
rT=reshape(rTrH(:,4),ni,nj);
ratio=reshape(rTrH(:,5),ni,nj);

% figure
% [c,h]=contour(T,Bp,rH,[0:10:90],'-k');
% clabel(c,h,[0:10:90],'fontsize',9)
% hold on
% [c,h]=contour(T,Bp,rT,[0:2:20],'-r');
% clabel(c,h,[0:2:20],'fontsize',9)
% set(gca,'XTick',[1 2 3 4 5 10 50 100])
% set(gca,'YTick',[50 100 200 300 400 500])
% set(gca,'XScale','log'); set(gca,'YScale','log')
% xlabel("planet's orbital period (day)",'fontsize',15)
% ylabel("planet's magnetic field (gauss)",'fontsize',15)
% title("Hill radius $r_H$ (black), magnetic truncation radius $r_T$ (red)",'fontsize',15,'interpreter','latex')
% print -dpdf rTrH.pdf

figure
%[c,h]=contourf(T,Bp,ratio,[0 1 2 5 10 15 20]);
[c,h]=contourf(T,Bp,ratio,[0:1:20]);
%clabel(c,h,[1 2 5 10 15 20],'fontsize',10)
set(h,'LineColor','none')
%colormap('jet')
colorbar
hold on
[c,h]=contour(T,Bp,ratio,[1:-1:0],'--r','LineWidth',3);
%hold on
%[c,h]=contour(T,Bp,ratio,[0:0.2:1],'-k');
%clabel(c,h,[0:0.2:1],'fontsize',10)
set(gca,'XTick',[1 2 3 4 5 10 50 100])
set(gca,'YTick',[50 100 200 300 400 500])
set(gca,'XScale','log'); set(gca,'YScale','log')
xlabel("planet's orbital period (day)",'fontsize',20)
ylabel("planet's magnetic field (gauss)",'fontsize',20)
title("ratio $r_H/r_T$",'fontsize',28,'interpreter','latex')
print -dpdf rTrH.pdf