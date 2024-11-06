% prova 2 MOM 
% project work 

clear all 
close all 
clc

%% Data 
kd=1/200; 
D=1.5; 
N=5000; 
xn=1000; 
x=1:1:5000; % chain length vector 

% initial distribution 
P0=zeros(1,N);
for n=1:N

    z=1/(D-1);
    y=1/D/xn*(z+1);       

    P0(n)= y^2/gamma(z)*n^(z-1)*exp(-y*n);
    
end

k=50; 

teta=0.01:50:6000;
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
cc=jet(6);
figure(1)
plot(teta,Pn_matrix(:,1).*1e4,"o",'Color',cc(1,:))
hold on
plot(teta,Pn_matrix(:,49).*1e4,"o",'Color',cc(2,:))
plot(teta,Pn_matrix(:,99).*1e4,"o",'Color',cc(3,:))
plot(teta,Pn_matrix(:,499).*1e4,"o",'Color',cc(4,:))
plot(teta,Pn_matrix(:,999).*1e4,"o",'Color',cc(5,:))
plot(teta,Pn_matrix(:,1999).*1e4,"o",'Color',cc(6,:))
xlabel('Dimensionless Time (Teta)')
ylabel('Normalized Concentration N*10^4')

