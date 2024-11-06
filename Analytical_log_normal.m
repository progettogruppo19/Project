% Log-normal distribution
% project work 

clear all 
close all 
clc

%% Data 
kd=1/200; 
D=1.5; 
N=5000; 
xn=1000; 
x=1:1:N; % chain length vector 

% initial distribution 
P0=zeros(1,N);
for nn=1:N
sigma=(log(D))^0.5;
mu=log(xn)-sigma^2/2;

P0(nn)=1/((2*pi)^0.5*nn*sigma)*exp(-((log(nn)-mu)^2)/(2*sigma^2));

end

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

%% plots

cc=jet(4);

plot(x,P0.*1e4,"o",'Color',cc(1,:))
hold on
plot(x,Pn_matrix(11,:).*1e4,"o",'Color',cc(2,:))
plot(x,Pn_matrix(26,:).*1e4,"o",'Color',cc(3,:))
plot(x,Pn_matrix(50,:).*1e4,"o",'Color',cc(4,:))

xlabel('Chain Length')
ylabel('Normalized Concentration N*10^4')
legend('initial distribution','teta= 100','teta=250','teta=500')