%pbe monodispersa+analitica monodispersa 5A

%% Introduction

clear
close
clc

%% Data

global N kd

kd = 1/200; %[1/s]
N = 1000;
chain_length = 1:1:1000;
D=1.5;
xn=50; 
x=1:1:N; % chain length vector 

% initial distribution 
P0=zeros(1,N);
P0(1,1000)=1;
k=50; 

teta=0.01:10:500;
tetaspan = 0:1:500; %[s]

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

[teta,P] = ode15s(@PBE, tetaspan, P0);
N=P./1;

%% plots
figure(1)
cc=jet(4);
%plot(x,P0.*1e4,"o",'Color',cc(1,:))

plot(x,Pn_matrix(6,:).*1e4,"o",'Color',cc(1,:))
hold on
plot(x,Pn_matrix(11,:).*1e4,"o",'Color',cc(2,:))
plot(x,Pn_matrix(26,:).*1e4,"o",'Color',cc(3,:))
plot(x,Pn_matrix(50,:).*1e4,"o",'Color',cc(4,:))

plot(chain_length, N(50,:).*1e4, 'LineWidth',1.6, 'Color','black')
hold on
plot(chain_length, N(100,:).*1e4, 'LineWidth',1.6)
plot(chain_length, N(250,:).*1e4, 'LineWidth',1.6)
plot(chain_length, N(500,:).*1e4, 'LineWidth',1.6)
axis([2.3 1000 0 750])
legend('Teta = 50', 'Teta = 100','Teta = 250','Teta = 500')
title('Monodisperse Distribution')
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