%% Introduction

clear
close
clc

options = odeset('AbsTol',1e-12, 'RelTol',1e-12);

%% Data

global kd N

kd = 1/200; %[1/s]
N = 200;
chain_length = 1:200;
xn = 50;

%Flory-Schultz distribution
 
P = 1 - 1/xn;
P0 = P.^(chain_length - 1).*(1-P);
 
tetaspan = 0:1:500; %[s]

%% Resolution

[teta,P] = ode15s(@PBE, tetaspan, P0, options);

%% Plots

figure(1)

plot(chain_length, P0.*1e4, 'LineWidth',1.6, 'Color','black')
hold on
plot(chain_length, P(10,:).*1e4, 'LineWidth',1.6)
plot(chain_length, P(25,:).*1e4, 'LineWidth',1.6)
plot(chain_length, P(45,:).*1e4, 'LineWidth',1.6)
axis([2.3 100 0 200])
legend('Teta = 0', 'Teta = 100','Teta = 250','Teta = 500')
title('Flory Distribution')
xlabel('Chain Length')
ylabel('Normalized Concentration N*10^4')

%% Function

function dPdteta = PBE(teta,P)

global N

%Initialisation

dPdteta = zeros(N,1);

%PBEs

for i = 3 : N

dPdteta(1) = sum(P(i)) + 2*P(2);

end

for n = 2 : N-1

dPdteta(n) = P(n+1) - P(n);

end

dPdteta(N) = -P(N);

end

