%% Introduction

clear
close
clc

%% Data 

global P0 maxCL x Pn_matrix teta_span

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

teta_span = [0 600];

%ODE

options = odeset('RelTol',1e-4, 'AbsTol',1e-4);
[teta,y] = ode15s(@MOM_fun,teta_span,initialcond, options);

t_ad = teta';
lambda0 = y(:,1);
lambda1 = y(:,2);
lambda2 = y(:,3);
M = y(:,4);

%% Plots 

figure(1)
plot(t_ad,lambda0',LineWidth=1.5)
title('Lambda 0')
xlabel('Dimensionless time (Teta)')
ylabel('Lambda 0')
xlim([0 600])

figure(2)
plot(t_ad,lambda1',LineWidth=1.5)
title('Lambda 1')
xlabel('Dimensionless time (Teta)')
ylabel('Lambda 1')
xlim([0 600])

figure(3)
plot(t_ad,lambda2',LineWidth=1.5)
title('Lambda 2')
xlabel('Dimensionless time (Teta)')
ylabel('Lambda 2')
xlim([0 600])

figure(4)
plot(t_ad,M',LineWidth=1.5)
xlabel('Dimensionless time (Teta)')
ylabel('Monomer concentration')
title('M')
xlim([0 600])

%% Functions 

function mom = MOM_fun(teta,y)

global P0 maxCL x teta_span Pn_matrix

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