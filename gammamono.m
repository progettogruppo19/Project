%gamma-mono, se aprite questo codice: va meglio nella prima parte ma penso
%scenda troppo verso la fine

%% Introduction

clear
close
clc

%% Data

global kd N

kd = 1/200; %[1/s]
chain_length = 1:1000; 
N = length(chain_length); 
D_gamma = 1.5;
D_mono=1;
xn = 1000;
 
z = 1/(D_gamma-1);
y = 1/D_gamma/xn*(z+1);
gamma = factorial(z - 1);
 
P0_gamma = y.^z./gamma.*chain_length.^(z-1).*exp(-y.*chain_length); 
P0_mono=zeros(1,N);
P0_mono(1,1000)=10;

P0=P0_gamma*0.94+P0_mono*0.06;
tetaspan = 0:1:500; %[s]

%% Resolution

[teta, P] = ode15s(@PBE, tetaspan, P0);

% Analytical solution

k = 50; 
teta_1 = 0.01:1:500;
Pn_total = zeros(length(teta_1), N+1);
for t = 1:length(teta_1)
    Pn_final = zeros(1, N+1); 
    for n = 2:N
        Pn_it = 0;
        for i = n:N
            Pn_old = Pn_it;
            if i < n + k + 1
                Pn_in = P0(i) * (teta_1(t))^(i-n) / factorial(i-n) * exp(-teta_1(t));
                Pn_it = Pn_in + Pn_old;
            else
                Pn_in = P0(i) / sqrt(2*pi*teta_1(t)) * exp(-(i-n-teta_1(t))^2 / (2*teta_1(t)));
                Pn_it = Pn_in + Pn_old;
            end
        end
        Pn_final(n) = Pn_it;
    end
    Pn_total(t, :) = Pn_final;
end

Pn_matrix = Pn_total(:, 2:end);

%% Plots
figure(1)
subplot(1,2,1);
plot(chain_length, P(51,:).*1e4, 'LineWidth',1.6, 'Color', 'black')
hold on
plot(chain_length,Pn_matrix(50,:).*1e4,"o", 'Color','black')
hold on
plot(chain_length, P(101,:) .* 1e4, 'LineWidth',1.2, 'Color','[0.3 0.2 0.8]')
plot(chain_length, Pn_matrix(100,:) .* 1e4, "o", 'Color','[0.3 0.2 0.8]')

plot(chain_length, P(251,:) .* 1e4, 'LineWidth',1.2, 'Color','[0.8 0.2 0.3]')
plot(chain_length, Pn_matrix(250,:) .* 1e4, "o", 'Color', '[0.8 0.2 0.3]')

plot(chain_length, P(501,:) .* 1e4, 'LineWidth',1.2, 'Color','[0.3 0.8 0.2]')
plot(chain_length, Pn_matrix(500,:) .* 1e4, "o", 'Color','[0.3 0.8 0.2]')

axis([50 1000 0 500]) 
legend('Teta = 50-Exact solution', 'Teta= 50-Analytical solution', 'Teta = 100-Exact Solution', 'Teta= 100-Analytical Solution', 'Teta = 250-Exact Solution', 'Teta = 250-Analytical Solution', 'Teta = 500 - Exact Solution', 'Teta = 500 - Analytical Solution')
title('Mixed model: Monodisperse+ Gamma Distribution')
xlabel('Chain Length')
ylabel('Normalized Concentration N*10^4')

subplot(1,2,2);
plot(chain_length, P(51,:).*1e4, 'LineWidth',1.6, 'Color', 'black')
hold on
plot(chain_length,Pn_matrix(50,:).*1e4,"o", 'Color','black')
hold on
plot(chain_length, P(101,:) .* 1e4, 'LineWidth',1.2, 'Color','[0.3 0.2 0.8]')
plot(chain_length, Pn_matrix(100,:) .* 1e4, "o", 'Color','[0.3 0.2 0.8]')

plot(chain_length, P(251,:) .* 1e4, 'LineWidth',1.2, 'Color','[0.8 0.2 0.3]')
plot(chain_length, Pn_matrix(250,:) .* 1e4, "o", 'Color', '[0.8 0.2 0.3]')

plot(chain_length, P(501,:) .* 1e4, 'LineWidth',1.2, 'Color','[0.3 0.8 0.2]')
plot(chain_length, Pn_matrix(500,:) .* 1e4, "o", 'Color','[0.3 0.8 0.2]')

axis([50 1000 0 50]) 
legend('Teta = 50-Exact solution', 'Teta= 50-Analytical solution', 'Teta = 100-Exact Solution', 'Teta= 100-Analytical Solution', 'Teta = 250-Exact Solution', 'Teta = 250-Analytical Solution', 'Teta = 500 - Exact Solution', 'Teta = 500 - Analytical Solution')
title('Mixed model: Monodisperse+ Gamma Distribution')
xlabel('Chain Length')
ylabel('Normalized Concentration N*10^4')

%% Function

function dPdteta = PBE(teta, P)

    global N

    % Inizializza dPdteta
    dPdteta = zeros(N,1);

    % PBEs
    for i = 3:N
        dPdteta(1) = sum(P(i)) + 2 * P(2);
    end

    for n = 2:N-1
        dPdteta(n) = P(n+1) - P(n);
    end

    dPdteta(N) = -P(N);

end
