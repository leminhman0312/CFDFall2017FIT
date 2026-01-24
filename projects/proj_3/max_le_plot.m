% cd Desktop/'Link to CFD_FALL_17'/HW3/

%{
filename = '005FTBS.dat';
FTBS005 = importdata(filename);
x = 1:1:71;
x = x.';
figure(1)
plot(x,FTBS005(:,1),'-xr')
hold on
box on
plot(x,FTBS005(:,2),'-sb')
hold on
plot(x,FTBS005(:,3),'-hg')
xlim([0 72])
xlabel('X [m]','fontsize',14);
ylabel('Displacement U [m]','fontsize',14);
title('$$\frac{\partial u}{\partial t} = {\alpha}\frac{\partial u}{\partial x}$$, FTBS explicit scheme,dt=0.005 sec','interpreter','latex','fontsize',18)
legend({'Initial condition', 'Half-way point','End condition'},'fontsize',14)





filename = '0025FTBS.dat';
FTBS0025 = importdata(filename);
x = 1:1:71;
x = x.';
figure(2)
plot(x,FTBS0025(:,1),'-xr')
hold on
box on
plot(x,FTBS0025(:,2),'-sb')
hold on
plot(x,FTBS0025(:,3),'-hg')
xlim([0 72])
xlabel('X [m]','fontsize',14);
ylabel('Displacement U [m]','fontsize',14);
title('$$\frac{\partial u}{\partial t} = {\alpha}\frac{\partial u}{\partial x}$$, FTBS explicit scheme,dt=0.0025 sec','interpreter','latex','fontsize',18)
legend({'Initial condition', 'Half-way point','End condition'},'fontsize',14)
%}


%{
filename = '0125FTBS.dat';
FTBS0125 = importdata(filename);
x = 1:1:71;
x = x.';
figure(3)
plot(x,FTBS0125(:,1),'-xr')
hold on
box on
plot(x,FTBS0125(:,2),'-sb')
hold on
plot(x,FTBS0125(:,3),'-hg')
xlim([0 72])
xlabel('X [m]','fontsize',14);
ylabel('Displacement U [m]','fontsize',14);
title('$$\frac{\partial u}{\partial t} = {\alpha}\frac{\partial u}{\partial x}$$, FTBS explicit scheme,dt=0.00125 sec','interpreter','latex','fontsize',18)
legend({'Initial condition', 'Half-way point','End condition'},'fontsize',14)

%}



%{
filename = '005Lax.dat';
Lax005 = importdata(filename);
x = 1:1:71;
x = x.';
figure(3)
plot(x,Lax005(:,1),'-xr')
hold on
box on
plot(x,Lax005(:,2),'-sb')
hold on
plot(x,Lax005(:,3),'-hg')
xlim([0 72])
xlabel('X [m]','fontsize',14);
ylabel('Displacement U [m]','fontsize',14);
title('$$\frac{\partial u}{\partial t} = {\alpha}\frac{\partial u}{\partial x}$$, Lax-Wendroff scheme,dt=0.005 sec','interpreter','latex','fontsize',18)
legend({'Initial condition', 'Half-way point','End condition'},'fontsize',14)




filename = '0025Lax.dat';
Lax0025 = importdata(filename);
x = 1:1:71;
x = x.';
figure(4)
plot(x,Lax0025(:,1),'-xr')
hold on
box on
plot(x,Lax0025(:,2),'-sb')
hold on
plot(x,Lax0025(:,3),'-hg')
xlim([0 72])
xlabel('X [m]','fontsize',14);
ylabel('Displacement U [m]','fontsize',14);
title('$$\frac{\partial u}{\partial t} = {\alpha}\frac{\partial u}{\partial x}$$, Lax-Wendroff scheme,dt=0.0025 sec','interpreter','latex','fontsize',18)
legend({'Initial condition', 'Half-way point','End condition'},'fontsize',14)




filename = '0125Lax.dat';
Lax0125 = importdata(filename);
x = 1:1:71;
x = x.';
figure(5)
plot(x,Lax0125(:,1),'-xr')
hold on
box on
plot(x,Lax0125(:,2),'-sb')
hold on
plot(x,Lax0125(:,3),'-hg')
xlim([0 72])
xlabel('X [m]','fontsize',14);
ylabel('Displacement U [m]','fontsize',14);
title('$$\frac{\partial u}{\partial t} = {\alpha}\frac{\partial u}{\partial x}$$, Lax-Wendroff scheme,dt=0.00125 sec','interpreter','latex','fontsize',18)
legend({'Initial condition', 'Half-way point','End condition'},'fontsize',14)
%}


%{
filename = '05implicit.dat';
implicit05 = importdata(filename);
x = 1:1:71;
x = x.';
figure(1)
plot(x,implicit05(:,1),'-xr')
hold on
box on
plot(x,implicit05(:,2),'-sb')
hold on
plot(x,implicit05(:,3),'-hg')
xlim([0 72])
xlabel('X [m]','fontsize',14);
ylabel('Displacement U [m]','fontsize',14);
title('$$\frac{\partial u}{\partial t} = {\alpha}\frac{\partial u}{\partial x}$$, BTCS implicit scheme,dt=0.005 sec','interpreter','latex','fontsize',18)
legend({'Initial condition', 'Half-way point','End condition'},'fontsize',14)



filename = '025implicit.dat';
implicit025 = importdata(filename);
x = 1:1:71;
x = x.';
figure(2)
plot(x,implicit025(:,1),'-xr')
hold on
box on
plot(x,implicit025(:,2),'-sb')
hold on
plot(x,implicit025(:,3),'-hg')
xlim([0 72])
xlabel('X [m]','fontsize',14);
ylabel('Displacement U [m]','fontsize',14);
title('$$\frac{\partial u}{\partial t} = {\alpha}\frac{\partial u}{\partial x}$$, BTCS implicit scheme,dt=0.0025 sec','interpreter','latex','fontsize',18)
legend({'Initial condition', 'Half-way point','End condition'},'fontsize',14)



filename = '0125implicit.dat';
implicit0125 = importdata(filename);
x = 1:1:71;
x = x.';
figure(3)
plot(x,implicit0125(:,1),'-xr')
hold on
box on
plot(x,implicit0125(:,2),'-sb')
hold on
plot(x,implicit0125(:,3),'-hg')
xlim([0 72])
xlabel('X [m]','fontsize',14);
ylabel('Displacement U [m]','fontsize',14);
title('$$\frac{\partial u}{\partial t} = {\alpha}\frac{\partial u}{\partial x}$$, BTCS implicit scheme,dt=0.00125 sec','interpreter','latex','fontsize',18)
legend({'Initial condition', 'Half-way point','End condition'},'fontsize',14)
%}





%{
filename = 'dt01burgers.dat';
u = importdata(filename);
x = 0:1:40;
x = x.';
u_table = u;
figure(7)
plot(x,u(:,1),'-xk','Linewidth',2)
hold on
box on
plot(x,u(:,2),'-^m')
hold on
plot(x,u(:,3),'-vb')
hold on 
plot(x,u(:,4),'->c')
hold on
plot(x,u(:,5),'-<g')
hold on
plot(x,u(:,6),'-sy')
hold on
plot(x,u(:,7),'-or')
xlabel('X [m]','fontsize',14);
ylabel('Displacement U [m]','fontsize',14);
title('$$\frac{\partial u}{\partial t} = -u\frac{\partial u}{\partial x}$$,MacCormack scheme,dt=0.1 sec','interpreter','latex','fontsize',18)
legend({'Initial condition', 't = 0.4 sec','t = 0.8 sec', 't = 1.2 sec', 't = 1.6 sec', 't = 2.0 sec', 't = 2.4 sec'},'fontsize',14)





filename = 'dt02burgers.dat';
u = importdata(filename);
x = 0:1:40;
x = x.';
figure(8)
plot(x,u(:,1),'-xk','Linewidth',2)
hold on
box on
plot(x,u(:,2),'-^m')
hold on
plot(x,u(:,3),'-vb')
hold on 
plot(x,u(:,4),'->c')
hold on
plot(x,u(:,5),'-<g')
hold on
plot(x,u(:,6),'-sy')
hold on
plot(x,u(:,7),'-or')
xlabel('X [m]','fontsize',14);
ylabel('Displacement U [m]','fontsize',14);
title('$$\frac{\partial u}{\partial t} = -u\frac{\partial u}{\partial x}$$,MacCormack scheme,dt=0.2 sec','interpreter','latex','fontsize',18)
legend({'Initial condition', 't = 0.4 sec','t = 0.8 sec', 't = 1.2 sec', 't = 1.6 sec', 't = 2.0 sec', 't = 2.4 sec'},'fontsize',14)
%}











