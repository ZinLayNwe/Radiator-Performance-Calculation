% Radiator Performance Calculation
% Written by DZLN & UPMT on 4-April-2024 ,  MATLAB version: R2023b
function [coolant_outT,T_Outlet_Air,eps] = evrad_calc(previous_value,TOutAir,Epsilon)
% Radiator Design Dimensions
Lrad = 0.63;                    % Length of Radiator, m
Hrad = 0.34;                    % Height of Radiator, m
Wrad = 0.016;                   % Width of Radiator, m

Htube = 0.002;                  % Height of coolant tube, m
Wtube = Wrad;                  % Width of coolant tube, m
Ntube = 35;                     % Number of coolant tube

Tfin = 0.00025;                 % Thickness of fin, m
Hfin = 0.00794;                 % Length or Heigth of fin, m
Wfin = Wtube;                   % Width of fin = width of tube
Fin_spacing = 0.002;            % fin - fin spacing per m

% Initial Condition
mdot_coolant = 0.094;           % Mass flow rate of coolant, kg/s
mdot_air = 0.2;                 % Mass flow rate of air, kg/s

T_inlet_coolant = 304.423;
T_inlet_air = 303.403;

T_Outlet_Coolant = previous_value;
T_Outlet_Air = TOutAir;

% Material Properties
% Coolant Properties
T_coolant = [293.15 313.15];
Rho_coolant = [1082.0 1069.0];
Cp_coolant = [3260.0 3340.0];
Mu_coolant = [0.00487 0.00257];
k_coolant = [0.402 0.398];
% Air Properties
T_Air = [300 350];
Rho_Air = [1.1614 0.995];
Cp_Air = [1007.0 1009.0];
Nu_Air = [1.589e-05 2.092e-5];
k_Air = [0.0263 0.03];
Pr_Air = [0.707 0.7];
% Aluminium Properties
T_Al = [200 400];
k_Al = [237 240];

% Functions Equations
% Outlet Temperatures 
Outlet_T_h = @(Thi,epsilon,Cmin,Ch,Tci) Thi-epsilon*Cmin/Ch*(Thi-Tci);
Outlet_T_c = @(Thi,epsilon,Cmin,Cc,Tci) Tci+epsilon*Cmin/Cc*(Thi-Tci);
% Log mean temperature difference
DTLM = @(x,y,z) ((x-z)-(x-y))/log((x-z)/(x-y));
% Mean Temperature
% Tmean = @(Tin,Tout) (Tin + Tout)/2;
% Specific Heat Capacity
% c = @(mdot,cp) mdot*cp;

% Calculation Start Here
Fin_density = 1/Fin_spacing;    % Fin density per meter
Nfin = Lrad * Fin_density;      % Number of fins per tude
Arad = Lrad * Hrad;             % Area of Radiator
% Coolant Tube side calculation
Atube = Wtube * Htube;          % cross sectional area of tude, m2 
Ptube = 2 * Wtube + 2 * Htube;  % Perimeter of tude, m
Dh = 4 * Atube / Ptube;         % Hydraulic diameter of tude, m
% Mean Temperature Calculation
T_mean_coolant = (T_inlet_coolant+T_Outlet_Coolant)/2;
T_mean_air = (T_inlet_air+T_Outlet_Air)/2;

Cp_coolant1 = interp1(T_coolant,Cp_coolant,T_mean_coolant);
Cp_Air1 = interp1(T_Air,Cp_Air,T_mean_air);
C_coolant = mdot_coolant * Cp_coolant1;
C_Air = mdot_air * Cp_Air1;
Cmax = max(C_Air,C_coolant);
Cmin = min(C_Air,C_coolant);
Cr = Cmin/Cmax; % Heat capacity ratio

Outlet_T_coolant = Outlet_T_h(T_inlet_coolant,Epsilon,Cmin,C_coolant,T_inlet_air);
Outlet_T_Air = Outlet_T_c(T_inlet_coolant,Epsilon,Cmin,C_Air,T_inlet_air);
% Mean Temperature by Heat capacity ratio
if Cr >= 0.5
    Thm = (T_inlet_coolant+Outlet_T_coolant)/2;
    Tcm = (T_inlet_air+Outlet_T_Air)/2;
elseif Cr < 0.5
    if sign(T_inlet_coolant-T_inlet_air)== 1
        if Cmax == C_coolant
            Thm = (T_inlet_coolant+Outlet_T_coolant)/2;
            DeltaTlm = DTLM(Thm,T_inlet_air,Outlet_T_Air);
            Tcm = Thm - DeltaTlm;
        end
        if Cmax == C_Air
            Tcm = (T_inlet_air+Outlet_T_Air)/2;
            DeltaTlm = DTLM(Tcm,T_inlet_coolant,Outlet_T_coolant);
            Thm = Tcm + DeltaTlm;
        end
    end
end
% Thermophysical Properties at Respective Mean Temperatures
% Coolant side
Q_coolant = mdot_coolant / interp1(T_coolant,Rho_coolant,Thm);
v_coolant = Q_coolant / (Ntube*Atube);                                       % velocity of coolant, m/s
rho_c_1 = interp1(T_coolant,Rho_coolant,Thm);
Mu_c_1 = interp1(T_coolant,Mu_coolant,Thm);
Re_coolant = v_coolant * rho_c_1 * Dh / Mu_c_1;
if Re_coolant <= 2000 && Wtube/Htube == 8
    Nu_coolant = 5.6;
end
K_C_1 = interp1(T_coolant,k_coolant,Thm);
h_coolant = Nu_coolant * K_C_1 / Dh;

% Air side
Q_Air = mdot_air / interp1(T_Air,Rho_Air,Tcm);
v_air = Q_Air / (Arad - (Ntube*Htube*Lrad));
Nu_A_1 = interp1(T_Air,Nu_Air,Tcm);
Re_air = v_air * Wfin / Nu_A_1;
Pr_A_1 = interp1(T_Air,Pr_Air,Tcm);
if Re_air <= 5e5 % Check for laminar flow
    Nu_Air_1 = 0.664 * sqrt(Re_air)* Pr_A_1^(1/3);
else % if turbulent flow
    % Add here equations for turbulent flow
end
K_A_1 = interp1(T_Air,k_Air,Tcm);
h_air = Nu_Air_1 * K_A_1 / Wtube;

% Fin Efficiency Calculations
k_Al = interp1(T_Al,k_Al,Thm);
m = sqrt(2*h_air/k_Al/Tfin);
Lc = Hfin + (Tfin/2);
etha_fin = tanh(m*Lc)/(m*Lc);
Afin_single = 2 * Wfin * Lc;
A_base_surface = 2 * Lrad * Wtube - Tfin * Wfin * Nfin;
A_Fin_base = Nfin * Afin_single + A_base_surface;
etha_overall = 1 - (Nfin*Afin_single/A_Fin_base)*(1-etha_fin);

A_external = A_Fin_base * Ntube;
A_internal = Ptube* Lrad*Ntube;

R_overall = (1/(etha_overall*h_air*A_external))+(1/(h_coolant*A_internal));	% Overall Thermal Resistance
UA = 1 / R_overall;                                               		% Overall Heat Transfer Coefficient
NTU = UA / Cmin;                                                          	 % Number of Transfer Units
eps = 1 - exp((1/Cr*NTU^0.22)*(exp(-Cr*NTU^0.78)-1));

qmax = Cmin * (T_inlet_coolant-T_inlet_air);
q_predicted = Epsilon * qmax;

coolant_outT = T_inlet_coolant - (q_predicted/C_coolant);
T_Outlet_Air = T_inlet_air + (q_predicted/C_Air);