%% Introduction

clear
close
clc

options = odeset('AbsTol',1e-10, 'RelTol',1e-10);

%% Data

global N

kd = 1/200; %[1/s]
N = 100;
chain_length = 1:100;
xn = 50;

%Poisson distribution

P0 = xn.^chain_length./factorial(chain_length).*exp(-xn);
 
tetaspan = 0:1:500; %[s]

%% Resolution

[teta,P] = ode15s(@PBE, tetaspan, P0, options);

%Analytical solution

k=50; 

teta_1=0.01:1:500;
Pn_total=zeros(length(teta_1),100+1);
for t=1:length(teta_1)
 Pn_final=zeros(1,100+1);
    for n=2:100 
 
    Pn_it=0;
    for i=n:100
    Pn_old=Pn_it;
        if i<n+k+1
        Pn_in=P0(i)*(teta_1(t))^(i-n)/(factorial(i-n))*exp(-teta_1(t));
        Pn_it=Pn_in+Pn_old;
        else
        Pn_in=P0(i)*1/(2*pi*teta_1(t))^0.5*exp(-(i-n-teta_1(t))^2/(2*teta_1(t)));
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

plot(chain_length, P0.*1e4, 'LineWidth',1.6, 'Color','black')
hold on
plot(chain_length, P(10,:).*1e4, 'LineWidth',1.6, 'Color','[0.3 0.2 0.8]')
plot(chain_length,Pn_matrix(10,:).*1e4,"o",  'Color','[0.3 0.2 0.8]')

plot(chain_length, P(25,:).*1e4, 'LineWidth',1.6, 'Color','[0.8 0.2 0.3]')
plot(chain_length,Pn_matrix(25,:).*1e4,"o", 'Color', '[0.8 0.2 0.3]')

plot(chain_length, P(45,:).*1e4, 'LineWidth',1.6, 'Color','[0.3 0.8 0.2]')
plot(chain_length,Pn_matrix(45,:).*1e4,"o", 'Color','[0.3 0.8 0.2]')

axis([3 100 0 600])
legend('Teta = 0', 'Teta = 100 - Exact Solution','Teta = 100 - Analytical Solution','Teta = 250 - Exact Solution','Teta = 250 - Analytical Solution','Teta = 500 - Exact Solution', 'Teta = 500 - Analytical Solution')
title('Poisson Distribution')
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