clear all 
clc
clf

filename1 = 'Yew_results_15000.dat';
U =importdata(filename1);

filename2 = 'Vee_results_15000.dat';
V =importdata(filename2);

filename3 = 'Pressure_results_15000.dat';
P =importdata(filename3);

sx = [0.0:0.05:0.5];
sy = [0.0:0.05:0.5];
[SX, SY] = meshgrid(sx,sy);

step_size = 0.1;
max_vertices = 2500;



dx = 0.00625;
nx = 81;
l = dx*(nx-1);

[x,y] = meshgrid(0:dx:0.5,0:dx:0.5);


% plotting
% figure(1)
streamline(x,y,U,V,SX,SY,[step_size,max_vertices]);
xlim([0 0.5])
ylim([0 0.5])
title('CAVITY FLOW STREAMLINE PLOT ','FontSize',15)
xlabel('X','FontSize',15);
ylabel('Y','FontSize',15)
xt = get(gca,'XTick');
set(gca, 'FontSize', 16);






%k = 2;
% figure(2)
%hq = quiver(x(1:k:end,1:k:end),y(1:k:end,1:k:end),U(1:k:end,1:k:end),V(1:k:end,1:k:end),3);
%xlim([0 0.5])
%ylim([0 0.5])
%title('CAVITY FLOW VECTOR PLOT (at 1/4 of the data points) ','FontSize',15)
%xlabel('X','FontSize',15);
%ylabel('Y','FontSize',15)
%xt = get(gca,'XTick');
%set(gca, 'FontSize', 16);






