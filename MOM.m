%% Introduction

clear
close
clc

%% Data 

global P0 maxCL x teta_span

kd = 1/200; 
D = 1.5;
xn = 1000; 
maxCL = 5000;
x = linspace(1,maxCL,maxCL);

%% Resolution

%Gamma distribution 

P0 = zeros(1,maxCL);

for n = 1:length(x)

    z = 1/(D-1);
    y = 1/D/xn*(z+1);
    gamma = factorial(z-1); 

    P0(n) = y^2/gamma*n^(z-1)*exp(-y*n);
   
end

%Inital conditions 

lambda0_in = sum(P0);
lambda1_in = sum(x.*P0);
lambda2_in = sum(x.^2.*P0);

initialcond =[lambda0_in lambda1_in lambda2_in 0];

teta_span = [0.01 600];

%ODE

options = odeset('RelTol',1e-4, 'AbsTol',1e-4);
[teta,y] = ode15s(@MOM_fun,teta_span,initialcond, options);

t_ad = teta';
lambda0 = y(:,1);
lambda1 = y(:,2);
lambda2 = y(:,3);
M = y(:,4);

yield_1 = 1 - lambda1./lambda1_in;
xn_1 = lambda1./lambda0;
xm_1 = lambda2./lambda1;
D_1 = xm_1./xn_1;

[teta_1,y_1] = ode15s(@MOM_fun_no_P2,teta_span,initialcond, options);

t_ad_1 = teta_1';
lambda0_1 = y_1(:,1);
lambda1_1 = y_1(:,2);
lambda2_1 = y_1(:,3);
M_1 = y(:,4);

yield_2 = 1 - lambda1_1/lambda1_in;
xn_2 = lambda1_1./lambda0_1;
xm_2 = lambda2_1./lambda1_1;
D_2 = xm_2./xn_2;

%% Plots 

figure(1)
plot(t_ad,lambda0',LineWidth=2)
hold on
plot(t_ad_1,lambda0_1',LineWidth=2)
title('Lambda 0')
xlabel('Dimensionless time (Teta)')
ylabel('Lambda 0')
axis([0 500 0 1.2])
legend('Solution considering P2', 'Solution without considering P2')

figure(2)
plot(t_ad,lambda1',LineWidth=2)
hold on
plot(t_ad_1,lambda1_1',"o",LineWidth=2)
title('Lambda 1')
xlabel('Dimensionless time (Teta)')
ylabel('Lambda 1')
xlim([0 500])
legend('Solution considering P2', 'Solution without considering P2')

figure(3)
plot(t_ad,lambda2',LineWidth=2)
hold on
plot(t_ad_1,lambda2_1',"o",LineWidth=2)
title('Lambda 2')
xlabel('Dimensionless time (Teta)')
ylabel('Lambda 2')
xlim([0 500])
legend('Solution considering P2', 'Solution without considering P2')

figure(4)
plot(t_ad,M',LineWidth=2)
xlabel('Dimensionless time (Teta)')
ylabel('Monomer concentration')
title('M')
xlim([0 500])

figure(5)
plot(t_ad, yield_1, LineWidth=2)
title('Yield with respect to Teta')
xlabel('Dimensionless time (Teta)')
ylabel('Yield')
axis([0 500 0 0.5])

figure(6)
plot(t_ad, xn_1, LineWidth=2)
hold on
plot(t_ad, xm_1, LineWidth=2)
title('Average Chain Length with respect to Teta')
xlabel('Dimensionless time (Teta)')
ylabel('Average Chain Length')
legend('Xn', 'Xm')
axis([0 500 600 1800])

figure(7)
plot(t_ad, D_1, LineWidth=2)
title('Dispersity with respect to Teta')
xlabel('Dimensionless time (Teta)')
ylabel('Dispersity')
axis([0 500 1.4 1.9])

%% Functions 

function mom = MOM_fun(teta,y)

global P0 maxCL x teta_span

lambda0 = y(1);
lambda1 = y(2);
lambda2 = y(3);
M = y(4);

k = 50;

Pn_final = zeros(1,maxCL+1);

for n = 2:maxCL 

    Pn_it = 0;  

    for ii = n:maxCL

        if ii <= k + n

            Pn_in = P0(ii) * (teta)^(ii - n) / factorial(ii - n) * exp(-teta);

        else

            Pn_in = P0(ii) / sqrt(2 * pi * teta) * exp(-(ii - n - teta)^2 / (2 * teta));

        end

        Pn_it = Pn_it + Pn_in;  
    end

    Pn_final(n) = Pn_it;  

end

Pn = Pn_final(2:end);

Pn_matrix = zeros(length(teta_span),length(Pn));

for i=1:length(teta)

 Pn_matrix(i,:)=Pn;

end

%Lambdas

lambda0 = sum(Pn);
lambda1 = sum(x.*Pn);
lambda2 = sum(x.^2.*Pn);

P2 = Pn(1);

%Equations 

dl0dteta = -P2;
dl1dteta = -lambda0-P2; 
dl2dteta = -2*lambda1+lambda0-P2; 
dMdteta = lambda0+P2; 

mom = [dl0dteta dl1dteta dl2dteta dMdteta]';

end

function mom = MOM_fun_no_P2(teta,y)

global P0 maxCL x teta_span

lambda0 = y(1);
lambda1 = y(2);
lambda2 = y(3);
M = y(4);

k = 50;

Pn_final = zeros(1,maxCL+1);

for n = 2:maxCL 

    Pn_it = 0;  

    for ii = n:maxCL

        if ii <= k + n

            Pn_in = P0(ii) * (teta)^(ii - n) / factorial(ii - n) * exp(-teta);

        else

            Pn_in = P0(ii) / sqrt(2 * pi * teta) * exp(-(ii - n - teta)^2 / (2 * teta));

        end

        Pn_it = Pn_it + Pn_in;  
    end

    Pn_final(n) = Pn_it;  

end

Pn = Pn_final(2:end);

Pn_matrix = zeros(length(teta_span),length(Pn));

for i=1:length(teta)

 Pn_matrix(i,:)=Pn;

end

%Lambdas

lambda0 = sum(Pn);
lambda1 = sum(x.*Pn);
lambda2 = sum(x.^2.*Pn);

%Equations 

dl0dteta = 0;
dl1dteta = -lambda0; 
dl2dteta = -2*lambda1+lambda0; 
dMdteta = lambda0; 

mom = [dl0dteta dl1dteta dl2dteta dMdteta]';

end