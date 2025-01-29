close all; clear; clc

load fort.3
x1=fort;
load fort.4
y1=fort;
load fort.31
x2=fort;
load fort.41
y2=fort;

figure
semilogy(x1(:,1),x1(:,2),'-k','linewidth',3)
hold on
semilogy(x1(:,1),x1(:,3),'-b','linewidth',3)
hold on
semilogy(x2(:,1),x2(:,3),'--b','linewidth',3)
hold on
semilogy(x1(:,1),x1(:,4),'-r','linewidth',3)
hold on
plot([3 3],[1d4 1d11],'--k','linewidth',1)
hold on
plot([40.0/24.0 40.0/24.0],[1d4 1d8],'--k','linewidth',1)
xlim([1 10])
ylim([1d4 1d11])
legend('$\Gamma_{\star p}^t$ with $\Omega_\star=3{\rm d}$','$\Gamma_{p\star}^t$ with $\Omega_p=40{\rm h}$','$\Gamma_{p\star}^t$ with $\Omega_p=3{\rm d}$','$\Gamma_{\star p}^m$ with $\Omega_\star=3{\rm d}$','corotation','fontsize',15,'location','southeast','interpreter','latex')
xlabel("planet's orbital period (day)", 'fontsize', 20)
ylabel("migration timescale (yr)", 'fontsize', 20)
title("planet's migration timescales", 'fontsize', 22)
print -dpdf migration1.pdf

figure
semilogy(y1(:,1),y1(:,2),'-k','linewidth',3)
hold on
semilogy(y1(:,1),y1(:,3),'-r','linewidth',3)
hold on
semilogy(y2(:,1),y2(:,2),'--k','linewidth',3)
hold on
semilogy(y2(:,1),y2(:,3),'--r','linewidth',3)
hold on
plot([8.7148 8.7148],[1d4 1d11],'--k','linewidth',1)
xlim([6 20])
ylim([1d8 1d12])
xticks([6 10 15 20 25 30])
xticklabels([6 10 15 20 25 30])
legend('$\Gamma_{pm}^t$ with $\Omega_p=40{\rm h}$','$\Gamma_{pm}^m$ with $\Omega_p=40{\rm h}$','$\Gamma_{pm}^t$ with $\Omega_p=3{\rm d}$','$\Gamma_{pm}^m$ with $\Omega_p=3{\rm d}$','fontsize',15,'location','southeast','interpreter','latex')
xlabel("moon's semi-major axis (R_J)", 'fontsize', 20)
ylabel("migration timescale (yr)", 'fontsize', 20)
title("moon's migration timescales", 'fontsize', 22)
print -dpdf migration2.pdf