clear;clc

cc= hsv(6); 
filepath= 'C:\Users\hxiang\OneDrive - UGent\Documents\PhD-data\phase field simulation\code\code_XH\PFM_Twinning\twinning simulation\6variants\PFM_twinning_6Variants_01\';
dataset= load([filepath 'VolumnFraction.mat']);
VolF= dataset.VolF;
    
f= figure('visible','on');
plot(VolF(:,1)* 1e-2, VolF(:,2), 'LineWidth',2, 'Color',cc(1,:)); hold on; plot(VolF(:,1)* 1e-2, VolF(:,3), 'LineWidth',2, 'Color',cc(2,:));hold on;
plot(VolF(:,1)* 1e-2, VolF(:,4), 'LineWidth',2, 'Color',cc(3,:)); hold on; plot(VolF(:,1)* 1e-2, VolF(:,5), 'LineWidth',2, 'Color',cc(4,:));hold on;
plot(VolF(:,1)* 1e-2, VolF(:,6), 'LineWidth',2, 'Color',cc(5,:)); hold on; plot(VolF(:,1)* 1e-2, VolF(:,7), 'LineWidth',2, 'Color',cc(6,:));grid on; hold off;
xlabel('Reduced simulation time (s)'); ylabel('Volume Fraction (%)');
legend({'V1', 'V2', 'V3', 'V4', 'V5', 'V6'},'FontSize',10,'NumColumns',1, 'TextColor','k'); legend boxoff; 
title(['Volume Fraction'],'Color','k');

filename= sprintf('VolumeFraction.tiff'); saveas(f,filename);f;clf