%% Introduction

clear
close
clc

%% Data

global N P_pbe kd

kd = 1/200; %[1/s]
chain_length = 1:5000; 
N = length(chain_length); 
D = 1.5;
xn = 1000;
 
z = 1/(D-1);
y = 1/D/xn*(z+1);
gamma_val = gamma(z); % Usa la funzione gamma per numeri non interi
 
P0 = y.^z./gamma_val.*chain_length.^(z-1).*exp(-y.*chain_length); 
lambda0_in = sum(P0);
lambda1_in = sum(chain_length.*P0);
lambda2_in = sum(chain_length.^2.*P0);

tetaspan = [0.01 600];

%% Risoluzione della distribuzione P usando PBE
[teta, P] = ode15s(@PBE, tetaspan, P0);
P_pbe=P(2);

%% Risoluzione delle variabili lambda usando lambda
initialcond = [lambda0_in lambda1_in lambda2_in 0];
[teta, y] = ode15s(@lambda, tetaspan, initialcond);

% Estrai i risultati
t_ad = teta';
lambda0 = y(:,1);
lambda1 = y(:,2);
lambda2 = y(:,3);
M = y(:,4);

% Grafico
figure(1)
plot(t_ad, lambda0)
xlabel('Dimensionless Time (teta)')
ylabel('Lambda0')

figure(2)
plot(t_ad, lambda1)
xlabel('Dimensionless Time (teta)')
ylabel('Lambda1')

figure(3)
plot(t_ad, lambda2)
xlabel('Dimensionless Time (teta)')
ylabel('Lambda2')

%% Funzioni

% Funzione per PBE
function dPdteta = PBE(teta, P)
    global N

    % Inizializza dPdteta
    dPdteta = zeros(N,1);

    % PBEs
    dPdteta(1) = sum(P(3:N)) + 2 * P(2);

    for n = 2:N-1
        dPdteta(n) = P(n+1) - P(n);
    end

    dPdteta(N) = -P(N);
end

% Funzione per lambda
function F = lambda(teta, y)
global P_pbe kd
    lambda0 = y(1);
    lambda1 = y(2);
    lambda2 = y(3);
    M = y(4);

    % Equazioni
    dl0dteta = -P_pbe;
    dl1dteta = -lambda0-P_pbe; 
    dl2dteta = -2*lambda1 +lambda0-P_pbe; 
    dMdteta = lambda0+P_pbe; 
    F = [dl0dteta; dl1dteta; dl2dteta; dMdteta];
end
