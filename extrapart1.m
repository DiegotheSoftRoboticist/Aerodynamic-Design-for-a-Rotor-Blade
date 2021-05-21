clear all
close all
clc

%Problem 5 
disp('Blade Element Momentum Iterator for a wind turbine')
disp('to be used to build up a whole blade BEM model')
disp('By Diego Ruiz & Ignacio Losada')


%Given: 

B = 3;
TSR = 8;
alpha = 6.1;
STA = 0;
Cl_max = 1.23;
R = 67.11;

%Solution. PROBLEM 4 Homework #5

r_R = (0.20:0.05:0.95);
LSR = TSR*r_R;

% if the twist angle is the sum of the pitch angle and the initial pitch
%angle where twist angle is zero, the initial pitch angle is -1.6.

alpha1=zeros(1,length(r_R));
alpha1(1,:) = 6.1;

%Solution. PROBLEM 1 Homework #6

phi_6 = (2/3)*atand(1./LSR);
r_6 = r_R*R;
c_6 = ((8*pi*r_6)/(B*Cl_max)).*(1-cosd(phi_6));
SP_6 = phi_6-alpha;
ST_6 = SP_6-SP_6(1,16)

%figure 6
Solidity = (B*c_6)./(2*pi*r_6);
a_6 =(4*sind(phi_6).^2);
a_6_b = (Solidity.*cosd(phi_6));
a_6 = 1./(1+(a_6)./a_6_b);
a(1:length(r_R)) = 0.3;

fig1 = figure;
c_6_R = c_6./R;
xlabel('r/R');ylabel('c/R');hold on;grid on
plot(r_R,c_6_R,'LineWidth',2);grid on
legend('with wake')

fig2 = figure;
grid on;hold on;
plot(r_R,ST_6,'LineWidth',2);grid on;
xlabel('r/R')
ylabel('Angle^\circ')
legend('Twist angle^\circ with wake')

fig3 = figure;
grid on;hold on;
plot(r_R,phi_6,'LineWidth',2);grid on;
xlabel('r/R')
ylabel('Angle^\circ')
legend('Angle of relative wind^\circ with wake');


fig4 = figure;
grid on;hold on;
plot(r_R,SP_6,'LineWidth',2);
xlabel('r/R')
ylabel('Angle^\circ')
legend('Section pitch angle^\circ with wake')



fig5 = figure;

grid on; hold on
plot(r_R,alpha1,'LineWidth',2);grid on
xlabel('r/R');
ylabel('angle of attack')
legend('Angle of Attack^\circ with wake')

fig6 = figure;


grid on;hold on;
plot(r_R,a_6,'LineWidth',2);grid on
legend('a with wake')
xlabel('r/R');
ylabel('Induction factor')


fig7 = figure;

a_prime(1:length(r_R)) = 0.0;
a_prime_6 = (1-3*a_6)./((4*a_6)-1);

grid on;hold on;
plot(r_R,a_prime_6,'LineWidth',2);grid on
xlabel('r/R');
ylabel('Induction factor')
legend('a_(_p_r_i_m_e_) with wake')


table = [r_6;c_6_R;ST_6;phi_6;SP_6;alpha1;a_6;a_prime_6]'

filename = 'extradesignoftheblade.xlsx';
xlswrite(filename,table)