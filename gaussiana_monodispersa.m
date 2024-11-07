%% Introduction

clear
close
clc

%% Data

global kd N

kd = 1/200; %[1/s]
N = 1000;
chain_length = 1:1000;

%Gamma distribution

D = 1.5;
xn = 1000;
z = 1/(D-1);
y = 1/D/xn*(z+1);
gamma = factorial(z - 1);
 
P0_gamma = y.^z./gamma.*chain_length.^(z-1).*exp(-y.*chain_length); 

%Monodisperse distribution

P0_mono = zeros(1,1000);
P0_mono(1000) = 1;

P0 = 0.06*P0_mono + 0.94*P0_gamma;

tspan = 0:10:500*200; %[s]

%% Resolution

[t,P] = ode15s(@PBE, tspan, P0);

teta = kd.*t;
N = P;

%Analytical solution

k=50; 

teta_1=0.01:5:500;
Pn_total=zeros(length(teta_1),1000+1);
for t=1:length(teta_1)
 Pn_final=zeros(1,1000+1);
    for n=2:1000 
 
    Pn_it=0;
    for i=n:1000
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

plot(chain_length, P(50,:).*1e4, 'LineWidth',1.6, 'Color','black')
hold on
plot(chain_length,Pn_matrix(6,:).*1e4,"o", 'Color','black')

plot(chain_length, P(100,:).*1e4, 'LineWidth',1.6, 'Color','[0.3 0.2 0.8]')
plot(chain_length,Pn_matrix(11,:).*1e4,"o", 'Color','[0.3 0.2 0.8]')

plot(chain_length, P(250,:).*1e4, 'LineWidth',1.6, 'Color','[0.8 0.2 0.3]')
plot(chain_length(1:3:end),Pn_matrix(26,1:3:end).*1e4,"o", 'Color','[0.8 0.2 0.3]')

plot(chain_length, P(500,:).*1e4, 'LineWidth',1.6, 'Color','[0.3 0.8 0.2]')
plot(chain_length(1:3:end),Pn_matrix(50,1:3:end).*1e4,"o", 'Color','[0.3 0.8 0.2]')
axis([5 1000 0 50])
legend('Teta = 50 - Exact Solution','Teta = 50 - Analytical Solution', 'Teta = 100 - Exact Solution','Teta = 100 - Analytical Solution','Teta = 250 - Exact Solution', 'Teta = 250 - Analytical Solution', 'Teta = 500 - Exact Solution', 'Teta = 500 - Analytical Solution')
title('Mono Dispersed Distribution')
xlabel('Chain Length')
ylabel('Normalized Concentration N*10^4')

%% Function

function dPdt = PBE(t,P)

global kd 

N = 1000;

%Initialisation

dPdt = zeros(N,1);

%PBEs

dPdt(1) = kd*sum(P(3:N)) + 2*kd*P(2);

for n = 2 : N-1

dPdt(n) = kd*(P(n+1) - P(n));

end

dPdt(N) = -kd*P(N);

end

