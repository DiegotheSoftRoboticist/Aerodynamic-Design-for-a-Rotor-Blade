%Problem 3.11


clear all
clc
close all


%Inicialitation. %Problem 6 Blade Shape Optimum Rotor with Wake Rotation.

% DATA for 1 N load.
filename ='designoftheblade.xlsx';
num1 = xlsread(filename);

r_R_311 = (0.05:0.1:0.95)';

B = 3; %# number of blades

TSR_311 = 9; % Tip Speed Ratio

LSR_311 = TSR_311*r_R_311;

% r_311 = num1(1:end,1);

R = 51;
r_311 = r_R_311*R;
Chord_311_R = num1(1:end,2);
Chord_311 = Chord_311_R.*R;
STA_311_deg = num1(1:end,3);
STA_311_rad=STA_311_deg*pi/180;
SP_311_deg = num1(1:end,5);
SP_311_rad = SP_311_deg*pi/180;
Cl_max_311 = 1.23;
Solidity_311 = B*Chord_311./(2*pi*r_311);


% if the twist angle is the sum of the pitch angle and the initial pitch
%angle where twist angle is zero, the initial pitch angle is -1.6

%Solution. PROBLEM 1 Homework #6

% STA_311_deg = [13;11;9;7;5;3.4;2.2;1.4;0.7;0.2]; % Section Twist Angle
% STA_311_rad=STA_311_deg*pi/180;
% SP_311_rad = STA_311_rad-((1.77*pi)/180);
% Chord_311 = [1;1;1;1;0.87;0.72;0.61;0.54;0.47;0.42];
% Solidity_311 = B*Chord_311./(2*pi*r_311);
% Cl_max_311 = 1;
% LSR_311 = TSR_311*r_R_311;





% c = 1; %chord
relax=0.1;
converge = 0.0001;
m=0;
ite = 1;
for TSR=1:17

    m=m+1;
    TSR_cp(m)=TSR;
    
    r_R = 0.05:0.1:0.95;
    LSR= (TSR_cp(m)*r_R)';
    
%INITITALIZATION.
%     a(i,1) = 1./(1+((4*sin(phi_6_rad).^2)./(Cl_max_6*Solidity_6.*cos(phi_6_rad))))
%     a_prime(i,1) = (1-3*a(i,1))./((4.*a(i,1))-1)
 
    error_axial = 1;
    error_rotation = 1;

%Iterations

    for i=1:10

        n = 1;

        a(1:10,1)=1/3;
        a(1:10,2)=0;
        a_prime(1:10,1)=0;
        a_prime(1:10,2)=0;
        while error_axial(n)>converge && error_rotation(n)>converge 
 
            
            
        phi_ite(i,ite) = atan(((1-a(i,ite))./(LSR(i,ite).*(1+a_prime(i,ite)))));

        F(i,ite) = (2/pi)*acos(exp(-(((B/2)*(1-r_R(ite,i)))'./(r_R(ite,i)'.*sin(phi_ite(i,ite))))));

        alpha(i,ite) = phi_ite(i,ite)-SP_311_rad(i,ite);

        alpha_deg(i,ite) = (alpha(i,ite)*180)/pi;

%         if alpha_deg(i,ite) >= -20 && alpha_deg(i,ite) <= 90
%               Cl_ite(i,ite) = -4.906789560853660E-10*(alpha_deg(i,ite)^6)+ 1.442768813298790E-07*(alpha_deg(i,ite)^5)...
%                   - 1.608493976865220E-05*(alpha_deg(i,ite)^4)+ 8.360188208672010E-04*(alpha_deg(i,ite)^3)...
%                   - 2.022769271535290E-02*(alpha_deg(i,ite)^2)+ 2.062537651900840E-01*(alpha_deg(i,ite))...
%                   + 0.105;
%               Cd_ite(i,ite) = +1.466647557343170E-10*(alpha_deg(i,ite)^6)- 3.605590522578590E-08*(alpha_deg(i,ite)^5)...
%                   + 3.234069226504270E-06*(alpha_deg(i,ite)^4)- 1.384701212572280E-04*(alpha_deg(i,ite)^3)...
%                   + 3.432617357852050E-03*(alpha_deg(i,ite)^2)- 2.068911716673940E-02*(alpha_deg(i,ite))...
%                   + 1.17e-02;
%             elseif alpha_deg(i,ite) > 90
%                   disp 'alpha > 90'
%             elseif alpha_deg(i,ite) >= -90 && alpha_deg(i,ite) < 20.2
%               Cl_ite(i,ite) = 0*(alpha_deg(i,ite)^6)- 5.1692355828E-09*(alpha_deg(i,ite)^5)...
%                   - 1.6355119057E-06*(alpha_deg(i,ite)^4)- 1.9133970760E-04*(alpha_deg(i,ite)^3)...
%                   - 9.5705220382E-03*(alpha_deg(i,ite)^2)- 1.8165835654E-01*(alpha_deg(i,ite))...
%                   - 1.6588452104;
%               Cd_ite(i,ite) = 0*(alpha_deg(i,ite)^6)- 1.1152917967E-08*(alpha_deg(i,ite)^5)...
%                   - 2.8195257909E-06*(alpha_deg(i,ite)^4)- 2.5518305432E-04*(alpha_deg(i,ite)^3)...
%                   - 1.0185092709E-02*(alpha_deg(i,ite)^2)- 2.1598008964E-01*(alpha_deg(i,ite))...
%                   - 1.6032899751;
%             end
if alpha_deg(i,ite) >= -90 && alpha_deg(i,ite) < -20.2
              Cl_ite(i,ite) = 0*(alpha_deg(i,ite)^6)- 5.1692355828E-09*(alpha_deg(i,ite)^5)...
                  - 1.6355119057E-06*(alpha_deg(i,ite)^4)- 1.9133970760E-04*(alpha_deg(i,ite)^3)...
                  - 9.5705220382E-03*(alpha_deg(i,ite)^2)- 1.8165835654E-01*(alpha_deg(i,ite))...
                  - 1.6588452104;
              Cd_ite(i,ite) = 0*(alpha_deg(i,ite)^6)- 1.1152917967E-08*(alpha_deg(i,ite)^5)...
                  - 2.8195257909E-06*(alpha_deg(i,ite)^4)- 2.5518305432E-04*(alpha_deg(i,ite)^3)...
                  - 1.0185092709E-02*(alpha_deg(i,ite)^2)- 2.1598008964E-01*(alpha_deg(i,ite))...
                  - 1.6032899751;
            elseif alpha_deg(i,ite) >= -20.2 && alpha_deg(i,ite) < 15.2 
              Cl_ite(i,ite) = + 5.8680602716E-08*(alpha_deg(i,ite)^6)+ 1.1307018778E-06*(alpha_deg(i,ite)^5)...
                  - 2.2279581308E-05*(alpha_deg(i,ite)^4)- 5.5537885629E-04*(alpha_deg(i,ite)^3)...
                  + 2.0710034350E-03*(alpha_deg(i,ite)^2)+ 1.2284570320E-01*(alpha_deg(i,ite))...
                  + 6.4951257976E-02;
              Cd_ite(i,ite) = - 1.2376869750E-08*(alpha_deg(i,ite)^6)- 1.0239769291E-07*(alpha_deg(i,ite)^5)...
                  + 6.4410496667E-06*(alpha_deg(i,ite)^4)+ 1.4537948353E-05*(alpha_deg(i,ite)^3)...
                  - 3.9178699176E-04*(alpha_deg(i,ite)^2)- 1.3992764977E-04*(alpha_deg(i,ite))...
                 + 1.4782368968E-02;
            elseif alpha_deg(i,ite) >= 15.2 && alpha_deg(i,ite) < 30.2
              Cl_ite(i,ite) = - 4.906789560853660E-10*(alpha_deg(i,ite)^6)+ 1.442768813298790E-07*(alpha_deg(i,ite)^5)...
                  - 1.608493976865220E-05*(alpha_deg(i,ite)^4)+ 8.360188208672010E-04*(alpha_deg(i,ite)^3)...
                  - 2.022769271535290E-02*(alpha_deg(i,ite)^2)+ 2.062537651900840E-01*(alpha_deg(i,ite))...
                  + 0.105;
              Cd_ite(i,ite) = + 1.466647557343170E-10*(alpha_deg(i,ite)^6)- 3.605590522578590E-08*(alpha_deg(i,ite)^5)...
                  + 3.234069226504270E-06*(alpha_deg(i,ite)^4)- 1.384701212572280E-04*(alpha_deg(i,ite)^3)...
                  + 3.432617357852050E-03*(alpha_deg(i,ite)^2)- 2.068911716673940E-02*(alpha_deg(i,ite))...
                  + 1.17e-02;
            elseif alpha_deg(i,ite) >= 30.2 && alpha_deg(i,ite) <= 90
              Cl_ite(i,ite) = 0*(alpha_deg(i,ite)^6)- 1.351495434687250E-08*(alpha_deg(i,ite)^5)...
                  + 4.323360471115850E-06*(alpha_deg(i,ite)^4)- 5.289720600966720E-04*(alpha_deg(i,ite)^3)...
                  + 3.018583104282480E-02*(alpha_deg(i,ite)^2)- 7.916138518339260E-01*(alpha_deg(i,ite))...
                  + 8.636443181225;
              Cd_ite(i,ite) = 0*(alpha_deg(i,ite)^6)+ 1.029729429071550E-08*(alpha_deg(i,ite)^5)...
                  - 2.480789504877290E-06*(alpha_deg(i,ite)^4)+ 2.042689103922690E-04*(alpha_deg(i,ite)^3)...
                  - 6.541889658905120E-03*(alpha_deg(i,ite)^2)+ 9.175984462513040E-02*(alpha_deg(i,ite))...
                  + 1.17e-02;
            elseif alpha_deg(i,ite) <-90 || alpha_deg(i,ite) > 90
                disp 'alpha = '; disp (alpha_deg(i,ite));
                disp ' i = '; disp (i);
            end

        Ct(i,ite) = Solidity_311(i,ite).*((1-a(i,ite)).^2).*(Cl_ite(i,ite).*cos(phi_ite(i,ite))+Cd_ite(i,ite).*sin(phi_ite(i,ite)))...
        ./(sin(phi_ite(i,ite)).^2);

    
            if Ct(i,ite)<0.96
        
                 a(i,ite+1) = 1./(1+((4.*F(i,ite).*sin(phi_ite(i,ite)).^2)./(Cl_ite(i,ite).*Solidity_311(i,ite).*cos(phi_ite(i,ite)))));
    
            elseif  Ct(i,ite)>0.96
        
                 a(i,ite+1) = (1./F(i,ite)).*(0.143+sqrt(0.0203-0.6427*(0.889-Ct(i,ite))));
            end
 
                 a_prime(i,ite+1) = 1./(((4.*F(i,ite).*cos(phi_ite(i,ite)))./(Cl_ite(i,ite).*Solidity_311(i,ite)))-1);

            n = n+1;
 
            error_axial(n) = abs((a(i,ite)-a(i,ite+1)));
            error_rotation(n) = abs((a_prime(i,ite)-a_prime(i,ite+1)));

            a_diff = a(i,ite+1)-(a(i,ite));
            a_prime_diff = ((a_prime(i,ite+1))-(a_prime(i,ite)));


a(i,ite) = a(i,ite)+relax*a_diff;
a_prime(i,ite) = a_prime(i,ite)+relax*a_prime_diff;






end

phi_cp(i,m)=phi_ite(i,ite);
cd_cp(i,m) = Cd_ite(i,ite);
cl_cp(i,m) = Cl_ite(i,ite);
LSR_cp(i,m) = LSR(i,ite);
F_cp(i,m) = F(i,ite);


end
end



Cp= (F_cp.*(sin(phi_cp).^2))...
.*(cos(phi_cp)-LSR_cp.*sin(phi_cp))...
.*(sin(phi_cp)+LSR_cp.*cos(phi_cp))...
.*((1-(cd_cp./cl_cp).*cot(phi_cp)).*LSR_cp.^2);


for m =1:length(TSR_cp)
  Cp_final(m) = (8/(TSR_cp(m)*10))*sum(Cp(:,m));
end 






plot(TSR_cp,Cp_final,'o','LineWidth',2);hold on
plot(TSR_cp,Cp_final)
title('Power Curve vs Tip Speed Ratio')
ylabel('Cp')
xlabel('Tip speed ratio');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%Problem #9 Sec%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%Initialization%%%%%
% 
nu = 0.8;   %drivetrain losses.
R = 51;     %rotor radio in meters.
rho = 1.225; %density at sea level in kg/m^3.



j=0;
i=0;
    
    
for U = 1:1:25
%for U = 2:2:16

    j=j+1;
  U_matrix(j)=U;
for  TSR=1:17
    i=i+1 
    
    Omega(i,j) = (TSR*U)/R;
    Omega_RPM(i,j)=(Omega(i,j)*60)/(2*pi);
    P(i,j) = Cp_final(1,i)*0.5*rho*(pi*R^2)*U^3;
 
end
i=0;
end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%g)

for i = 1:length(P(1,:))
    
    P_max(i) = max(P(:,i))
    n = find(P(:,i) == P_max(i))
    Omega_RPM_max(i)=Omega_RPM(n,i);
     
    Omega_max(i)=Omega(n,i);

end
TSR_max = (Omega_max*R)/U

%  U = (0.5:0.1:1)*11;
%  %U = 2:2:16
% table_g = [U;P_max;TSR_max;Omega_max]'
% T = table(table_g)
% writetable(T,'Table g).xlsx') 

%Part h)

% n = 0;
% for Omega_h = [50.93,61.85,72.76]
%     n = n+1;
%     for i = 1:length(P(1,:))
%      P_H(n,i) = interp1(Omega_RPM(:,i),P(:,i),Omega_h)   
%     Omega_RPM_H(i,n)=Omega_h;
%     end
% end
% 
% for j=1:length(U)
% for i =1:3
% TSR_H(i,j) = (((2*pi)/60)*Omega_RPM_H(1,i)*R)/U(1,j)
% end
% end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig2 = figure

hold on
 for j= 1:length(Omega_RPM(1,:))
 
    plot(Omega_RPM(:,j),P(:,j),'LineWidth',1.2);
    plot(Omega_RPM_max(j),P_max(j),'ro','LineWidth',2);
 end
  hold off
  
title('Power vs. Angular Velocity')
ylabel('Power [Watts]')
xlabel('\omega [RPM]');

legend('2 m/s','Maximum for 2 m/s','4 m/s','Maximum for 4 m/s','6 m/s'...
    ,'Maximum for 6 m/s','8 m/s','Maximum for 8 m/s','10 m/s',...
    'Maximum for 8 m/s','12 m/s','Maximum for 12 m/s','14 m/s',...
    'Maximum for 14 m/s','16 m/s','Maximum for 16 m/s');

% %%%%%%%%%%%%%%%%%%%%%PART H plot
% fig4 = figure
% 
% hold on
%  for j= 1:length(Omega_RPM(1,:))
%  
%     plot(Omega_RPM(:,j),P(:,j),'LineWidth',0.9);
%     plot(Omega_RPM_H(j,:),P_H(:,j),'*--','LineWidth',1.5);
%  end
%   hold off
%   
% title('Power vs. Angular Velocity')
% ylabel('Power [Watts]')
% xlabel('\omega [RPM]');
% 
% legend('2 m/s','\omega_1 = \omega_o_p_t , \omega_2 = 0.85\omega_o_p_t , \omega_3 = 0.7\omega_o_p_t'...
%     ,'4 m/s','\omega_1 = \omega_o_p_t , \omega_2 = 0.85\omega_o_p_t , \omega_3 = 0.7\omega_o_p_t','6 m/s'...
%     ,'\omega_1 = \omega_o_p_t , \omega_2 = 0.85\omega_o_p_t , \omega_3 = 0.7\omega_o_p_t','8 m/s',...
%     '\omega_1 = \omega_o_p_t , \omega_2 = 0.85\omega_o_p_t , \omega_3 = 0.7\omega_o_p_t','10 m/s',...
%     '\omega_1 = \omega_o_p_t , \omega_2 = 0.85\omega_o_p_t , \omega_3 = 0.7\omega_o_p_t','12 m/s',...
%     '\omega_1 = \omega_o_p_t , \omega_2 = 0.85\omega_o_p_t , \omega_3 = 0.7\omega_o_p_t','14 m/s',...
%     '\omega_1 = \omega_o_p_t , \omega_2 = 0.85\omega_o_p_t , \omega_3 = 0.7\omega_o_p_t','16 m/s',...
%     '\omega_1 = \omega_o_p_t , \omega_2 = 0.85\omega_o_p_t , \omega_3 = 0.7\omega_o_p_t');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%f)
% P_opt = 0.9*194266; %Watts
% Omega_RPM_des = 72.76; %RPM
% U_des = 11; %m/s
% 
% 
% P_m = P/P_opt;
% Omega_m = Omega_RPM/Omega_RPM_des;
% U_m = U/U_des;
% 
% fig3 = figure
% 
%   hold on
%  for j= 1:length(Omega_RPM(1,:))
%      
%     plot(Omega_m(:,j),P_m(:,j),'LineWidth',1.2);
%     
%  end
%   hold off
%   
% title('SCIG')
% ylabel('P_m')
% xlabel('\omega_m');
% legend('v_w = 0.5','v_w = 0.6','v_w = 0.7','v_w =0.8','v_w = 0.9','v_w = 1','v_w =1.27','v_w = 1.46');
% hold on
% plot([1 1],[0 1],'r','LineWidth',2.2)

