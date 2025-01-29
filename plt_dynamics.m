clear; close all; clc

% initial am=6R_J
load dyn1.dat % planet's orbital period=4 d
load dyn2.dat % planet's orbital period=10 d
load dyn3.dat % planet's orbital period=100 d
figure
semilogx(dyn1(:,1),dyn1(:,3),'-k',dyn1(:,1),dyn1(:,4),'--k',dyn2(:,1),dyn2(:,3),'-r',dyn2(:,1),dyn2(:,4),'--r',dyn3(:,1),dyn3(:,3),'-b',dyn3(:,1),dyn3(:,4),'--b','linewidth',2)
hold on
plot([10 1e8],[0.1196093642E+01 0.1196093642E+01],':k','linewidth',1)
xlim([0 1e7])
ylim([0 20])
yticks([0:2:20])
text(20, 17, 'black: 4 days \newline {\color{red}red}: 10 days \newline {\color{blue}blue}: 100 days', 'fontsize', 18)
xlabel('time (yr)','fontsize',22)
ylabel('distance from planet (R_J)','fontsize',20)
title('r_{cor} (solid) and a_m (dashed)','fontsize',20)
print -dpng dyn1.png

% initial am=10R_J
load dyn4.dat % planet's orbital period=4 d
load dyn5.dat % planet's orbital period=10 d
load dyn6.dat % planet's orbital period=100 d
figure
semilogx(dyn4(:,1),dyn4(:,3),'-k',dyn4(:,1),dyn4(:,4),'--k',dyn5(:,1),dyn5(:,3),'-r',dyn5(:,1),dyn5(:,4),'--r',dyn6(:,1),dyn6(:,3),'-b',dyn6(:,1),dyn6(:,4),'--b','linewidth',2)
hold on
plot([10 1e8],[0.1196093642E+01 0.1196093642E+01],':k','linewidth',1)
xlim([0 1e8])
ylim([0 22])
yticks([0:2:22])
text(20, 19, 'black: 4 days \newline {\color{red}red}: 10 days \newline {\color{blue}blue}: 100 days', 'fontsize', 18)
xlabel('time (yr)','fontsize',22)
ylabel('distance from planet (R_J)','fontsize',20)
title('r_{cor} (solid) and a_m (dashed)','fontsize',20)
print -dpng dyn2.png