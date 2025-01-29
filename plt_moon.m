clear; clc; close all

for i=1:100
map(i,1)=0.7; map(i,2)=0.7; map(i,3)=0.7;
end

load moon2.dat %Mm=1.d26
ni=200; nj=200;
Pp=reshape(moon2(:,1),nj,ni+1);
Pm=reshape(moon2(:,2),nj,ni+1);
tau2=reshape(moon2(:,4),nj,ni+1);
beta2=reshape(moon2(:,5),nj,ni+1);

figure
hold 
[c,h]=contourf(Pp,Pm,beta2,[-10:0.5:5]);
set(h,'LineColor','none')
colorbar
[c,h]=contour(Pp,Pm,beta2,[0 5],'-r','LineWidth',3);
[c,h]=contour(Pp,Pm,tau2,[1:1:1],'--r','LineWidth',3);
%plot([0.186140E+03,0.186140E+03],[2,0.186140E+03],':k','linewidth',2)
xlim([2 200])
ylim([2 200])
set(gca,'XTick',[2 3 4 5 10 50 100 150 200])
set(gca,'YTick',[2 3 4 5 10 50 100 150 200])
set(gca,'XScale','log'); set(gca,'YScale','log')
text(140,3.5,'A','FontSize',20)
text(15,4,'B','FontSize',20)
text(60,18,'C','FontSize',20)
text(8,50,'D','FontSize',20)
xlabel("planet's orbital period (day)",'fontsize',20)
ylabel("moon's orbital period (day)",'fontsize',20)
title("$\log(\tau_{\rm cor}/\tau_{\rm mig})$",'fontsize',28,'interpreter','latex')
print -dpdf tau1.pdf

load moon1.dat %Mm=1.d25
load moon3.dat %Mm=1.d27
tau1=reshape(moon1(:,4),nj,ni+1);
beta1=reshape(moon1(:,5),nj,ni+1);
tau3=reshape(moon3(:,4),nj,ni+1);
beta3=reshape(moon3(:,5),nj,ni+1);
figure
hold
[c,h]=contour(Pp,Pm,beta1,[0 5],'-k','LineWidth',2); % Mm=1.d25
[c,h]=contour(Pp,Pm,tau1,[1:1:1],'--k','LineWidth',2);
[c,h]=contour(Pp,Pm,beta2,[0 5],'-r','LineWidth',2); % Mm=1.d26
[c,h]=contour(Pp,Pm,tau2,[1:1:1],'--r','LineWidth',2);
[c,h]=contour(Pp,Pm,beta3,[0 5],'-b','LineWidth',2); % Mm=1.d27
[c,h]=contour(Pp,Pm,tau3,[1:1:1],'--b','LineWidth',2);
%plot([0.186140E+03,0.186140E+03],[2,0.186140E+03],':k','linewidth',2)
plot([2,200],[2,200],':k','LineWidth',1)
xlim([2 200])
ylim([2 200])
set(gca,'XTick',[2 3 4 5 10 50 100 150 200])
set(gca,'YTick',[2 3 4 5 10 50 100 150 200])
set(gca,'XScale','log'); set(gca,'YScale','log')
text(3,80,'black: M_m=10^{25}g \newline {\color{red}red}: M_m=10^{26}g \newline {\color{blue}blue}: M_m=10^{27}g','fontsize',18)
xlabel("planet's orbital period (day)",'fontsize',20)
ylabel("moon's orbital period (day)",'fontsize',20)
title("$\tau_{\rm cor}/\tau_{\rm mig}=1$",'fontsize',28,'interpreter','latex')
print -dpdf tau2.pdf