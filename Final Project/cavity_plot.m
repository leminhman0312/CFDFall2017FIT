clear all 
clc

filename1 = 'Yew_results.dat';
U =importdata(filename1);

filename2 = 'Vee_results.dat';
V =importdata(filename2);

filename3 = 'Pressure_results.dat';
P =importdata(filename3);

sx = [0.0:0.05:0.5];
sy = [0.0:0.05:0.5];
[SX, SY] = meshgrid(sx,sy);

step_size = 0.1;
max_vertices = 2500;



dx = 0.00625;
nx = 11;
l = dx*(nx-1);

x = 0:dx:l;
y = 0:dx:l;


%plotting

% streamline(x,y,U.',V.',SX,SY,[step_size,max_vertices]);

quiver(x,y,U.',V.');




