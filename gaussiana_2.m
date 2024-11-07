%% Introduction

clear
close
clc

%% Data

global kd N

kd = 1/200; %[1/s]
N = 5000;
chain_length = 1:5000;

%Gamma distribution 1

D = 1.2;
xn_1 = 500;
z = 1/(D-1);
y_1 = 1/D/xn_1*(z+1);
gamma = factorial(5 - 1);
 
P0_1 = y_1.^z./gamma.*chain_length.^(z-1).*exp(-y_1.*chain_length); 

%Gamma distribution 2

xn_2 = 2000;
y_2 = 1/D/xn_2*(z+1);
P0_2 = y_2.^z./gamma.*chain_length.^(z-1).*exp(-y_2.*chain_length); 

P0 = 0.2*P0_1 + 0.8*P0_2;

tspan = 0:10:500*200; %[s]

%% Resolution

[t,P] = ode15s(@PBE, tspan, P0);

teta = kd.*t;
N = P;

%Analytical solution

k=50; 

teta_1=0.01:5:500;
Pn_total=zeros(length(teta_1),5000+1);
for t=1:length(teta_1)
 Pn_final=zeros(1,5000+1);
    for n=2:5000 
 
    Pn_it=0;
    for i=n:5000
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

plot(chain_length, P0 .* 1e4, 'LineWidth',1.6, 'Color','black')
hold on
plot(chain_length, P(10,:) .* 1e4, 'LineWidth',1.6, 'Color','[0.3 0.2 0.8]')
plot(chain_length(1:10:end), Pn_matrix(11,1:10:end) .* 1e4, "o", 'Color','[0.3 0.2 0.8]')

plot(chain_length, P(25,:) .* 1e4, 'LineWidth',1.6, 'Color','[0.8 0.2 0.3]')
plot(chain_length(1:10:end), Pn_matrix(26,1:10:end) .* 1e4, "o", 'Color', '[0.8 0.2 0.3]')

plot(chain_length, P(45,:) .* 1e4, 'LineWidth',1.6, 'Color','[0.3 0.8 0.2]')
plot(chain_length(1:7:end), Pn_matrix(50,1:7:end) .* 1e4, "o", 'Color','[0.3 0.8 0.2]')

axis([50 5000 0 5]) 
legend('Teta = 0', 'Teta = 100 - Exact Solution', 'Teta = 100 - Analytical Solution', 'Teta = 250 - Exact Solution', 'Teta = 250 - Analytical Solution', 'Teta = 500 - Exact Solution', 'Teta = 500 - Analytical Solution')
title('4 Gamma Distributions')
xlabel('Chain Length')
ylabel('Normalized Concentration N*10^4')

%% Function

function dPdt = PBE(t,P)

global kd N

%Initialisation

dPdt = zeros(N,1);

%PBEs

dPdt(1) = kd*sum(P(3:N)) + 2*kd*P(2);

for n = 2 : N-1

dPdt(n) = kd*(P(n+1) - P(n));

end

dPdt(N) = -kd*P(N);

end

