%% Introduction

clear
close
clc

%% Data

global N kd

kd = 1/200; %[1/s]
N = 1000;
chain_length = 1:1:1000;

%Monodisperse distribution

P0 = zeros(1,999);
P0 = [P0 1];

tetaspan = 0:1:500; %[s]

%% Resolution

[teta,P] = ode15s(@PBE, tetaspan, P0);

%% Plots

figure(1)

plot(chain_length, P0.*1e4, 'LineWidth',1.6, 'Color','black')
hold on
plot(chain_length, P(100,:).*1e4, 'LineWidth',1.6)
plot(chain_length, P(250,:).*1e4, 'LineWidth',1.6)
plot(chain_length, P(500,:).*1e4, 'LineWidth',1.6)
axis([4 1000 0 100])
legend('Teta = 0', 'Teta = 100','Teta = 250','Teta = 500')
title('Log-normal Distribution')
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