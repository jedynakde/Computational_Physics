% test 1 for T >1

T =    [2.04,2.01,1.99,1.97,2   ,2.02 ,1.98 ,2.07  ,2     ,2];
P =    [2.1 ,3   ,4.5 ,6.6 ,9.47,15.16,22.83,33.57 ,49.366,121.423];
Pnid = [1.05,1.8 ,3.2 ,5.1 ,7.66,13.14,20.65,31.075,46.77 ,118.42];
rho =  [0.5 ,0.6 ,0.7 ,0.8 ,0.9 ,1.0  ,1.1  ,1.2   ,1.3   ,1.4];
% T variation = 0.2
% P variation = 0.1
% Pnid variation = 0.1

plot(rho,T,rho,P,rho,Pnid);
%itle("Isotherm of T ~ 2 from rho 0.5 to 1.4");
xlabel ("Density(rho)");
legend ("Temperature","Pressure","Ideal Pressure");


% test 2 

T =    [0.05 ,0.05   ,0.05 ,0.05   ,0.05   ,0.05   ,0.05  ,0.05   ,0.05    ,0.05    ,0.05];
P =    [0.031,-0.0049,-0.01,-.00042,0.0018   ,-0.0035,0.0026,-0.0023,-0.00189,0.001186,2.15*10^-5];
Pnid = [0.005,-0.023 ,-0.025,-0.01 ,-0.0059,-0.0088,-0.0016,-0.0054,-0.0022 ,0.00017,-0.000485];
rho =  [0.5 ,0.4     ,0.3   ,0.2   ,0.15   ,0.1    ,0.08   ,0.06   ,.04   ,0.02     ,   0.01];
% T variation = 0.2
% P variation = 0.1
% Pnid variation = 0.1
figure;
plot(rho,T,rho,P,rho,Pnid);
title("Isotherm of T ~ 0.05 from rho 0.5 to 0.01");
xlabel ("Density(rho)");
legend ("Temperature","Pressure","Ideal Pressure");





 