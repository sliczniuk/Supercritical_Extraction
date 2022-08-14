clc, close all
clear all
%addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.2');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

%% Parameters
nstages = 5;                                                   %

C0fluid = 0;                                            % Extractor initial concentration of extract
% Fluid phase kg / m^3
C0solid = 32;                                           % Solid phase kg / m^3
V       = 0.01/4;                                       % Volume of the extractor  [m3] - Vargas
epsi    = 0.70;                                         % Porosity [-] - Vargas
dp      = 0.001;                                        % Diameter of the particle [m] - Vargas
L       = 0.70;                                         % Length of the extractor [m] - Vargas
rho_s   = 1087.2;                                       % Densisty of the solid phase [kg / m^3] - Vargas
km      = 0.29;                                         % Partition coefficient (?)
mi      = 1/5;                                          % Geometric shape coefficient (?) - Vargas
Ti      = [313.15, 323.15, 333.15, 343.15];
Di      = [ 3.03E-11, 3.66E-11, 1.97E-10, 2.73E-10];
T_kp    = [313.15, 323.15, 333.15, 343.15];
dkp     = [8.13E-01 6.94E-01 1.04E-01 6.67E-02];


Tc      = 304.1;                                        % Critical temperature [K]
Pc      = 73.8;                                         % Critical pressure [bar]
omega   = 0.228;                                        % Acentric factor [-]
R       = 8.314e-5;                                     % Universal gas constant, [m3-bar/K-mol]

kappa   = 0.37464 + 1.54226 * omega - 0.26992 * omega^2;

MW      = 44e-3;                                        % Molar weight of CO2 [Kg/mol]

[EA_Di, betah_Di] = Arrhenius_Di(Ti,Di,R);
[EA_km, betah_km] = Arrhenius_km(T_kp,dkp,R,rho_s);

CP_0 =   4.18686;
CP_A =  19.80;        %4.5980;
CP_B =   7.344E-2;    %0.0125;
CP_C = - 5.602E-5;    %2.86E-06;
CP_D =   1.715E-8;    %-2.70E-09;

cpSolid = 1.5E3;      % J / K / Kg

% axial diffusion correlation
a_axial = 0.1;
b_axial = 0.011;
c_axial = 0.48;

% Heat Conductivity
A1_cond = -105.161;
A2_cond =  0.9007;
A3_cond =  0.0007;
A4_cond =  3.50E-15;
A5_cond =  3.76E-10;
A6_cond =  0.7500;
A7_cond =  0.0017;

% Viscosity
A1_visc = -1.146067E-1;
A2_visc =  6.978380E-7;
A3_visc =  3.976765E-10;
A4_visc =  6.336120E-2;
A5_visc = -1.166119E-2;
A6_visc =  7.142596E-4;
A7_visc =  6.519333E-6;
A8_visc = -3.567559E-1;
A9_visc =  3.180473E-2;

%                 1        2    3    4    5  6     7   8   9   10  11  12   13   14   15      16       17    18    19    20    21    22      23        24      25        26       27      28       29        30        31       32       33       34      35       36        37      38       39       40      41        42       43   
parameters = {nstages, C0solid, V, epsi, dp, L, rho_s, km, mi, Tc, Pc, R, kappa, MW, EA_Di, betah_Di, CP_0, CP_A, CP_B, CP_C, CP_D, EA_km, betah_km, cpSolid, a_axial, b_axial, c_axial, A1_cond, A2_cond, A3_cond, A4_cond, A5_cond, A6_cond, A7_cond, A1_visc, A2_visc, A3_visc, A4_visc, A5_visc, A6_visc, A7_visc, A8_visc, A9_visc};


%% dummy parameters
clc

T_check = 0.8*Tc:0.05:1.1*Tc;
P_check = 0.8*Pc:0.05:1.2*Pc;


Z   = nan(numel(T_check), numel(P_check), 3 );
RHO = nan(numel(T_check), numel(P_check), 3 );
CP  = nan(numel(T_check), numel(P_check), 3 );
MU  = nan(numel(T_check), numel(P_check), 3 );
KT  = nan(numel(T_check), numel(P_check), 3 );


for i = 1:numel(T_check)
    for j = 1:numel(P_check)

        T = T_check(i);
        P = P_check(j);
        
        Tr      = T ./ Tc;
        a       = 0.4572350 .* R.^2 .* Tc.^2 ./ Pc;
        b       = 0.0777961 .* R    .* Tc    ./ Pc;

        alpha = (1 + kappa .* (1 - sqrt(Tr))).^2;

        A = a .* alpha .* P ./ R.^2 ./ T.^2;
        B = b .* P ./ R ./ T;

        z = roots([1, - (1 - B), (A - 2.*B - 3.*B.^2), - ( A .* B - B.^2 - B.^3) ]);

        z = z(real(z)>0&imag(z)==0);

        if numel(z) == 1
            Z(i,j,1) = z;
        else
            Z(i,j,2) = max(z);
            Z(i,j,3) = min(z); 
        end
        
        rho = full(rhoPB_Comp(T, P, z, parameters));
        if numel(z) == 1
            RHO(i,j,1) = rho;
        else
            RHO(i,j,2) = max(rho);
            RHO(i,j,3) = min(rho);
        end

        if numel(z) == 1
            CP(i,j,1) = SpecificHeatComp(T, P, z, rho, parameters);
        else
            CP(i,j,2) = SpecificHeatComp(T, P, max(z), max(rho), parameters);
            CP(i,j,3) = SpecificHeatComp(T, P, min(z), min(rho), parameters);
        end

        mu = Viscosity(T,P,parameters);
        if numel(z) == 1
            MU(i,j,1) = mu;
        else
            MU(i,j,2) = max(mu);
            MU(i,j,3) = min(mu);
        end

        if numel(z) == 1
            KT(i,j,1) = HeatConductivity_Comp(T, P, z, rho, parameters);
        else
            KT(i,j,2) = HeatConductivity_Comp(T, P, max(z), max(rho), parameters);
            KT(i,j,3) = HeatConductivity_Comp(T, P, min(z), min(rho), parameters);
        end
        
    end
end


%% Plotting

% Compressibility

figure(1)
set(gcf,'PaperOrientation','landscape')
subplot(1,3,1); contourf(P_check,T_check,squeeze(Z(:,:,1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
hcb=colorbar;
title(hcb,'Z');
title('Compressibility at one-phase region',' ')
xlabel('Pressure [bar]')
ylabel('Temperature [K]')

subplot(1,3,2); contourf(P_check,T_check,squeeze(Z(:,:,2)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
hcb=colorbar;
title(hcb,'Z');
title('Compressibility at two-phase region:',' the biggest root')
xlabel('Pressure [bar]')
ylabel('Temperature [K]')

subplot(1,3,3); contourf(P_check,T_check,squeeze(Z(:,:,3)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
hcb=colorbar;
title(hcb,'Z');
title('Compressibility at two-phase region:',' the smallest root')
xlabel('Pressure [bar]')
ylabel('Temperature [K]')

print(gcf, '-dpdf', '-fillpage', 'Compressibility.pdf'); close all

% Rho

figure(2)
set(gcf,'PaperOrientation','landscape')
subplot(1,3,1); contourf(P_check,T_check,squeeze(RHO(:,:,1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
hcb=colorbar;
title(hcb,'\rho_{CO_2}');
title('Density at one-phase region',' ')
xlabel('Pressure [bar]')
ylabel('Temperature [K]')

subplot(1,3,2); contourf(P_check,T_check,squeeze(RHO(:,:,2)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
hcb=colorbar;
title(hcb,'\rho_{CO_2}');
title('Density at two-phase region:',' the biggest root')
xlabel('Pressure [bar]')
ylabel('Temperature [K]')

subplot(1,3,3); contourf(P_check,T_check,squeeze(RHO(:,:,3)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
hcb=colorbar;
title(hcb,'\rho_{CO_2}');
title('Density at two-phase region:',' the smallest root')
xlabel('Pressure [bar]')
ylabel('Temperature [K]')

print(gcf, '-dpdf', '-fillpage', 'RHO.pdf'); close all

% CP

figure(3)
mCP = log10(CP);
set(gcf,'PaperOrientation','landscape')
subplot(1,3,1); contourf(P_check,T_check,squeeze(mCP(:,:,1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
hcb=colorbar;
title(hcb,'log_{10}(CP)');
title('Specific heat at one-phase region:',' ')
xlabel('Pressure [bar]')
ylabel('Temperature [K]')

subplot(1,3,2); contourf(P_check,T_check,squeeze(mCP(:,:,2)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
hcb=colorbar;
title(hcb,'log_{10}(CP)');
title('Specific heat at two-phase region:',' the biggest root')
xlabel('Pressure [bar]')
ylabel('Temperature [K]')

subplot(1,3,3); contourf(P_check,T_check,squeeze(mCP(:,:,3)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
hcb=colorbar;
title(hcb,'log_{10}(CP)');
title('Specific heat at two-phase region:',' the smallest root')
xlabel('Pressure [bar]')
ylabel('Temperature [K]')

print(gcf, '-dpdf', '-fillpage', 'CP.pdf'); close all

% mu

figure(4)
set(gcf,'PaperOrientation','landscape')
subplot(1,3,1); contourf(P_check,T_check,squeeze(MU(:,:,1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
hcb=colorbar;
title(hcb,'\mu');
title('Viscosity at one-phase region:',' ')
xlabel('Pressure [bar]')
ylabel('Temperature [K]')

subplot(1,3,2); contourf(P_check,T_check,squeeze(MU(:,:,2)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
hcb=colorbar;
title(hcb,'\mu');
title('Viscosity at two-phase region:',' the biggest root')
xlabel('Pressure [bar]')
ylabel('Temperature [K]')

subplot(1,3,3); contourf(P_check,T_check,squeeze(MU(:,:,3)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
hcb=colorbar;
title(hcb,'\mu');
title('Viscosity at two-phase region:',' the smallest root')
xlabel('Pressure [bar]')
ylabel('Temperature [K]')

print(gcf, '-dpdf', '-fillpage', 'MU.pdf'); close all


% kt

figure(5)
set(gcf,'PaperOrientation','landscape')
subplot(1,3,1); contourf(P_check,T_check,squeeze(KT(:,:,1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
hcb=colorbar;
title(hcb,'kt');
title('Thermal conductivity at one-phase region:',' ')
xlabel('Pressure [bar]')
ylabel('Temperature [K]')

subplot(1,3,2); contourf(P_check,T_check,squeeze(KT(:,:,2)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
hcb=colorbar;
title(hcb,'kt');
title('Thermal conductivity at two-phase region:',' the biggest root')
xlabel('Pressure [bar]')
ylabel('Temperature [K]')

subplot(1,3,3); contourf(P_check,T_check,squeeze(KT(:,:,3)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
hcb=colorbar;
title(hcb,'kt');
title('Thermal conductivity at two-phase region:',' the smallest root')
xlabel('Pressure [bar]')
ylabel('Temperature [K]')

print(gcf, '-dpdf', '-fillpage', 'KT.pdf'); close all
































