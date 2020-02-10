clear all 
clc
x = 0.2*(0:30);
y = 0.2*(0:20);


%Point Gauss Seidel

h1 = figure(1);
filename = 'ResultsPointGS.dat';
delimiterIn = '\t';
PointGaussSeidel = importdata(filename, delimiterIn);
contour(x,y,PointGaussSeidel)
h = colorbar;
set(h,'Ydir','reverse');
title('Solving $$\bigg[\frac{\partial^2 \Psi}{\partial x^2} + \frac{\partial^2 \Psi}{\partial y^2}\bigg]$$ = 0 using Point Gauss Seidel','interpreter','latex','fontsize',18);
xlabel('X [feet]','fontsize',14')
ylabel('Y [feet]','fontsize',14')

%{
%Line Gauss Seidel
h2 = figure(2);
filename = 'ResultsLineGS.dat';
delimiterIn = '\t';
LineGS = importdata(filename, delimiterIn);
contour(x,y,LineGS);
h = colorbar;
set(h,'Ydir','reverse');
title('Solving $$\bigg[\frac{\partial^2 \Psi}{\partial x^2} + \frac{\partial^2 \Psi}{\partial y^2}\bigg]$$ = 0 using Line Gauss Seidel','interpreter','latex','fontsize',18);
xlabel('X [feet]','fontsize',14)
ylabel('Y [feet]','fontsize',14)
%}


%{
%Point SOR
h3 = figure(3);
filename = 'ResultsPSOR.dat';
delimiterIn = '\t';
PointSOR = importdata(filename, delimiterIn);
contour(x,y,PointSOR);
h = colorbar;
set(h,'Ydir','reverse');
title('Solving $$\bigg[\frac{\partial^2 \Psi}{\partial x^2} + \frac{\partial^2 \Psi}{\partial y^2}\bigg]$$ = 0 using Point SOR','interpreter','latex','fontsize',18);
xlabel('X [feet],'fontsize',14')
ylabel('Y [feet],'fontsize',14')
%}




h4 = figure(4);
filename = 'testPSOR.dat';
delimiterIn = '\t';
TestPointSOR = importdata(filename, delimiterIn);
wRangePointSOR = TestPointSOR(1:20,1);
iterRangePointSOR = TestPointSOR(1:20,2);
plot(iterRangePointSOR,wRangePointSOR);
xlim([0 12000]);
title('Iteration vs relaxation parameter for Point SOR')
xlabel('Iteration','fontsize',14);
ylabel('{\omega}','fontsize',18);
set(get(gca,'ylabel'),'rotation',360);






h5 = figure(5);
filename = 'testLSOR.dat';
delimiterIn = '\t';
TestLineSOR = importdata(filename, delimiterIn);
wRangeLineSOR = TestLineSOR(1:14,1);
iterRangeLineSOR = TestLineSOR(1:14,2);
plot(iterRangeLineSOR,wRangeLineSOR);
%xlim([0 12000]);
title('Iteration vs relaxation parameter for Line SOR')
xlabel('Iteration','fontsize',14);
ylabel('{\omega}','fontsize',18);
set(get(gca,'ylabel'),'rotation',360);






%{
% This is for the PSOR with good value of W 
% Run one then change values 
%Point SOR
h6 = figure(6);
filename = 'omegaPSOR.dat';
delimiterIn = '\t';
PointSOR = importdata(filename, delimiterIn);
contour(x,y,PointSOR);
h = colorbar;
set(h,'Ydir','reverse');
title('Solving $$\bigg[\frac{\partial^2 \Psi}{\partial x^2} + \frac{\partial^2 \Psi}{\partial y^2}\bigg]$$ = 0 using Point SOR at w = 1.7','interpreter','latex','fontsize',18);
xlabel('X [feet]','fontsize',14)
ylabel('Y [feet]','fontsize',14)
%}





%{
% This is for the LSOR with good value of W 
% Run one then change values 
%Point SOR
%h7 = figure(7);
filename = 'omegaLSOR.dat';
delimiterIn = '\t';
LineSOR = importdata(filename, delimiterIn);
contour(x,y,LineSOR);
h = colorbar;
set(h,'Ydir','reverse');
title('Solving $$\bigg[\frac{\partial^2 \Psi}{\partial x^2} + \frac{\partial^2 \Psi}{\partial y^2}\bigg]$$ = 0 using Line SOR at w = 1.7','interpreter','latex','fontsize',18);
xlabel('X [feet]','fontsize',14)
ylabel('Y [feet]','fontsize',14)

%}














