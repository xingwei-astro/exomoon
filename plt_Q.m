clear; clc; close all
load Q.dat

figure
plot(Q(:,1),Q(:,6),'--r','linewidth',2)
hold on
plot(Q(:,1),Q(:,7),'--b','linewidth',2)

figure
semilogy(Q(:,1),Q(:,3))