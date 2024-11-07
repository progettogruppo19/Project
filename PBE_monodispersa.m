%% Introduction

clear
close
clc

%% Data

global N

kd = 1/200; %[1/s]
N = 1000;
chain_length = 1:1000;
x=1:0:N; % chain length vector 

%Monodisperse distribution

P0 = zeros(1,1000);
P0(1000) = 1;

tetaspan = 0:1:500; %[s]

%% Resolution

[teta,P] = ode15s(@PBE, tetaspan, P0);

%Analytical solution

k=50; 

teta=0.01:10:500;
Pn_total=zeros(length(teta),N+1);
for t=1:length(teta)
 Pn_final=zeros(1,N+1);
    for n=2:N 
 
    Pn_it=0;
    for i=n:N
    Pn_old=Pn_it;
        if i<n+k+1
        Pn_in=P0(i)*(teta(t))^(i-n)/(factorial(i-n))*exp(-teta(t));
        Pn_it=Pn_in+Pn_old;
        else
        Pn_in=P0(i)*1/(2*pi*teta(t))^0.5*exp(-(i-n-teta(t))^2/(2*teta(t)));
        Pn_it=Pn_in+Pn_old;
        end 
    end
 Pn_final(n)=Pn_it;
end
Pn_total(t,:)=Pn_final;
end

Pn_matrix=Pn_total(:,2:end);

%% Plots

figure(1)

plot(chain_length, P(50,:).*1e4, 'LineWidth',1.6, 'Color','black')
hold on
plot(chain_length,Pn_matrix(6,:).*1e4,"o", 'Color','black')

plot(chain_length, P(100,:).*1e4, 'LineWidth',1.6, 'Color','[0.3 0.2 0.8]')
plot(chain_length,Pn_matrix(11,:).*1e4,"o", 'Color','[0.3 0.2 0.8]')

plot(chain_length, P(250,:).*1e4, 'LineWidth',1.6, 'Color','[0.8 0.2 0.3]')
plot(chain_length,Pn_matrix(26,:).*1e4,"o", 'Color','[0.8 0.2 0.3]')

plot(chain_length, P(500,:).*1e4, 'LineWidth',1.6, 'Color','[0.3 0.8 0.2]')
plot(chain_length(1:2:end),Pn_matrix(50,1:2:end).*1e4,"o", 'Color','[0.3 0.8 0.2]')
axis([5 1000 0 750])
legend('Teta = 50 - Exact Solution','Teta = 50 - Analytical Solution', 'Teta = 100 - Exact Solution','Teta = 100 - Analytical Solution','Teta = 250 - Exact Solution', 'Teta = 250 - Analytical Solution', 'Teta = 500 - Exact Solution', 'Teta = 500 - Analytical Solution')
title('Mono Dispersed Distribution')
xlabel('Chain Length')
ylabel('Normalized Concentration N*10^4')

%% Function

function dPdteta = PBE(teta,P)

global N

%Initialisation

dPdteta = zeros(N,1);

%PBEs

dPdteta(1) = sum(P(3:N)) + 2*P(2);

for n = 2 : N-1

dPdteta(n) = P(n+1) - P(n);

end

dPdteta(N) = -P(N);

end