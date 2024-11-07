%% Introduction

clear
close
clc

%% Data

global kd N

kd = 1/200; %[1/s]
N = 5000;
chain_length = 1:5000;

%Gamma distribution

D = 1.5;
xn = 1000;
z = 1/(D-1);
y = 1/D/xn*(z+1);
gamma = factorial(z - 1);
 
P0 = y.^z./gamma.*chain_length.^(z-1).*exp(-y.*chain_length); 

tspan = 0:10:1200000; %[s]

%% Resolution

[t,P] = ode15s(@PBE, tspan, P0);

teta = kd.*t;
N = P;

%Analytical solution

k=50; 

teta_1=0.01:50:6000;
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

plot(chain_length(2), P0(2)*1e4, 'o', 'LineWidth',1.6)
hold on

title('Chain Length Distribution at Teta = 0')
xlabel('Chain Length')
ylabel('Normalized Concentration N*10^4')

plot(chain_length(50), P0(50)*1e4, 'o', 'LineWidth',1.6)

plot(chain_length(100), P0(100)*1e4, 'o', 'LineWidth',1.6)

plot(chain_length(500), P0(500)*1e4, 'o', 'LineWidth',1.6)

plot(chain_length(1000), P0(1000)*1e4, 'o', 'LineWidth',1.6)

plot(chain_length(2000), P0(2000)*1e4, 'o', 'LineWidth',1.6)

hold on
plot(chain_length, P0.*1e4, 'LineWidth',1.6, 'Color','black')
legend('P 2', 'P 50', 'P 100', 'P 500', 'P 1000', 'P 2000')
hold off

cc = jet(6);

figure(2)

plot(teta,N(:,2).*1e4, 'LineWidth',1.6, 'Color',cc(1,:))
hold on
plot(teta_1,Pn_matrix(:,1).*1e4,"o",'Color',cc(1,:))
plot(teta,N(:,50).*1e4, 'LineWidth',1.6, 'Color',cc(2,:))
plot(teta_1,Pn_matrix(:,49).*1e4,"o",'Color',cc(2,:))
plot(teta,N(:,100).*1e4, 'LineWidth',1.6, 'Color',cc(3,:))
plot(teta_1,Pn_matrix(:,99).*1e4,"o",'Color',cc(3,:))
plot(teta,N(:,500).*1e4, 'LineWidth',1.6, 'Color',cc(4,:))
plot(teta_1,Pn_matrix(:,499).*1e4,"o",'Color',cc(4,:))
plot(teta,N(:,1000).*1e4, 'LineWidth',1.6, 'Color',cc(5,:))
plot(teta_1,Pn_matrix(:,999).*1e4,"o",'Color',cc(5,:))
plot(teta,N(:,2000).*1e4, 'LineWidth',1.6, 'Color',cc(6,:))
plot(teta_1,Pn_matrix(:,1999).*1e4,"o",'Color',cc(6,:))
axis([0 6000 0 8])
xlabel('Dimensionless Time (Teta)')
ylabel('Normalized Concentration N*10^4')
legend('P0 2 Exact', 'P0 2 Analytical', 'P0 50 Exact', 'P0 50 Analytical','P0 100 Exact', 'P0 100 Analytical','P0 500 Exact', 'P0 500 Analytical','P0 1000 Exact','P0 1000 Analytical','P0 2000 Exact', 'P0 2000 Analytical')

%% Function

function dPdt = PBE(t,P)

global kd N

%Initialisation

dPdt = zeros(N,1);

%PBEs

for i = 3 : N

dPdt(1) = kd*sum(P(i)) + 2*kd*P(2);

end

for n = 2 : N-1

dPdt(n) = kd*(P(n+1) - P(n));

end

dPdt(N) = -kd*P(N);

end

