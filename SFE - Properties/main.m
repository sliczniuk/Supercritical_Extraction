clc, close all
clear all
addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.5');
%addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
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


Tc      = 304.2;                                        % Critical temperature [K]
Pc      = 73.765;                                       % Critical pressure [bar]
omega   = 0.225;                                        % Acentric factor [-]
R       = 83.14;                                        % Universal gas constant, [cm3-bar/mol-K]

kappa   = 0.37464 + 1.5422 * omega - 0.2699 * omega^2;

MW      = 44.01;                                        % Molar weight of CO2 [Kg/mol]

[EA_Di, betah_Di] = Arrhenius_Di(Ti,Di,R);
[EA_km, betah_km] = Arrhenius_km(T_kp,dkp,R,rho_s);

CP_0 =  4.1868;
CP_A =  4.7280;        %4.5980;
CP_B =  0.0175;    %0.0125;
CP_C = -1.34E-05;    %2.86E-06;
CP_D =  4.10E-09;    %-2.70E-09;

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

correlation_viscosity    = {'Amooey', 'Fenghour', 'Laesecke'};
correlation_conductivity = {'Amooey', 'Bahadori', 'Jarrahian', 'Rostami', 'Rostamian', 'Huber'};

T_check = linspace(Tc,1.1*Tc,100);
P_check = linspace(Pc,1.5*Pc,100);
%T_check = 50+273;
%P_check = 90;


Z   = nan(numel(P_check), numel(T_check), 3 );
RHO = nan(numel(P_check), numel(T_check), 3 );
CP  = nan(numel(P_check), numel(T_check), 3 );
MU  = nan(numel(P_check), numel(T_check), 3, numel(correlation_viscosity) );
KT  = nan(numel(P_check), numel(T_check), 3, numel(correlation_conductivity) );

for i = 1:numel(P_check)
    if mod(i,round(numel(P_check)/100)) == 0
        clc
        fprintf('%f %%\n',i/numel(P_check)*100)
    end

    for j = 1:numel(T_check)

        T = T_check(j);
        P = P_check(i);

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
            Z(i,j,1) = max(z);
            Z(i,j,3) = min(z);
        end

        rho = full(rhoPB_Comp(T, P, z, parameters)) ;

        if numel(z) == 1
            RHO(i,j,1) = rho;
        else
            RHO(i,j,1) = max(rho);
            RHO(i,j,3) = min(rho);
        end

        if numel(z) == 1
            CP(i,j,1) = SpecificHeatComp(T, P, z, rho, parameters);
        else
            CP(i,j,1) = SpecificHeatComp(T, P, max(z), max(rho), parameters);
            CP(i,j,3) = SpecificHeatComp(T, P, min(z), min(rho), parameters);
        end

        for k=1:numel(correlation_viscosity)
            correlation_mass  = correlation_viscosity{k};
            if numel(z) == 1
                MU(i,j,1,k) = Viscosity(T,P,rho,parameters,correlation_mass);
            else
                MU(i,j,2,k) = Viscosity(T,P,max(rho),parameters,correlation_mass);
                MU(i,j,3,k) = Viscosity(T,P,min(rho),parameters,correlation_mass);
            end
        end

        for k=1:numel(correlation_conductivity)
            correlation_heat = correlation_conductivity{k};
            if numel(z) == 1
                KT(i,j,1, k) = HeatConductivity_Comp(T, P, z, rho, parameters, correlation_heat);
            else
                KT(i,j,2, k) = HeatConductivity_Comp(T, P, max(z), max(rho), parameters, correlation_heat);
                KT(i,j,3, k) = HeatConductivity_Comp(T, P, min(z), min(rho), parameters, correlation_heat);
            end
        end

    end
end

%% Plotting

%{

% Compressibility

figure(1);
set(gcf,'PaperOrientation','landscape', 'visible','off')
subplot(1,3,1); contourf(T_check,P_check,squeeze(Z(:,:,1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
xline(Tc,':','color','w')
yline(Pc,':','color','w')
hcb=colorbar;
title(hcb,'Z');
title('Compressibility 1-phase region',' ')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,2); contourf(T_check,P_check,squeeze(Z(:,:,2)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
xline(Tc,':','color','w')
yline(Pc,':','color','w')
hcb=colorbar;
title(hcb,'Z');
title('Compressibility 2-phase region:','the biggest root')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,3); contourf(T_check,P_check,squeeze(Z(:,:,3)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
xline(Tc,':','color','w')
yline(Pc,':','color','w')
hcb=colorbar;
title(hcb,'Z');
title('Compressibility 2-phase region:','the smallest root')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

print(gcf, '-dpdf', '-fillpage', 'Compressibility.pdf'); close all

% Rho

figure(2)
set(gcf,'PaperOrientation','landscape', 'visible','off')
subplot(1,3,1); contourf(T_check,P_check,squeeze(RHO(:,:,1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
xline(Tc,':','color','w')
yline(Pc,':','color','w')
hcb=colorbar;
title(hcb,'\rho_{CO_2} [kg/m^3]');
title('Density 1-phase region',' ')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,2); contourf(T_check,P_check,squeeze(RHO(:,:,2)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
xline(Tc,':','color','w')
yline(Pc,':','color','w')
hcb=colorbar;
title(hcb,'\rho_{CO_2} [kg/m^3]');
title('Density 2-phase region:','the biggest root')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,3); contourf(T_check,P_check,squeeze(RHO(:,:,3)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
xline(Tc,':','color','w')
yline(Pc,':','color','w')
hcb=colorbar;
title(hcb,'\rho_{CO_2} [kg/m^3]');
title('Density 2-phase region:','the smallest root')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

print(gcf, '-dpdf', '-fillpage', 'RHO.pdf'); close all

% CP

figure(3)
mCP = log10(CP);
set(gcf,'PaperOrientation','landscape', 'visible','off')
subplot(1,3,1); contourf(T_check,P_check,squeeze(mCP(:,:,1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
xline(Tc,':','color','w')
yline(Pc,':','color','w')
hcb=colorbar;
title(hcb,'log_{10}(CP)');
title('Specific heat 1-phase region:',' ')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,2); contourf(T_check,P_check,squeeze(mCP(:,:,2)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
xline(Tc,':','color','w')
yline(Pc,':','color','w')
hcb=colorbar;
title(hcb,'log_{10}(CP)');
title('Specific heat 2-phase region:','the biggest root')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,3); contourf(T_check,P_check,squeeze(mCP(:,:,3)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
xline(Tc,':','color','w')
yline(Pc,':','color','w')
hcb=colorbar;
title(hcb,'log_{10}(CP)');
title('Specific heat 2-phase region:','the smallest root')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

print(gcf, '-dpdf', '-fillpage', 'CP.pdf'); close all

% mu

figure(4)
set(gcf, 'visible','off')

for i=0:numel(correlation_viscosity)-1
    subplot(numel(correlation_viscosity),3,3*i+1); contourf(T_check,P_check,squeeze(MU(:,:,1,i+1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
    xline(Tc,':','color','w')
    yline(Pc,':','color','w')
    hcb=colorbar;
    title(hcb,'\mu[Pa \times s]');
    title('Viscosity 1-phase region:',' ')
    
    ylabel({correlation_viscosity{i+1};'Pressure [bar]'})
    xlabel('Temperature [K]')
    
    subplot(numel(correlation_viscosity),3,3*i+2); contourf(T_check,P_check,squeeze(MU(:,:,2,i+1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
    xline(Tc,':','color','w')
    yline(Pc,':','color','w')
    hcb=colorbar;
    title(hcb,'\mu[Pa \times s]');
    title('Viscosity 2-phase region:','the biggest root')
    ylabel('Pressure [bar]')
    xlabel('Temperature [K]')
    
    subplot(numel(correlation_viscosity),3,3*i+3); contourf(T_check,P_check,squeeze(MU(:,:,3,i+1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
    xline(Tc,':','color','w')
    yline(Pc,':','color','w')
    hcb=colorbar;
    title(hcb,'\mu[Pa \times s]');
    title('Viscosity 2-phase region:','the smallest root')
    ylabel('Pressure [bar]')
    xlabel('Temperature [K]')
end
print(gcf, '-dpdf', '-fillpage', 'MU.pdf'); close all


% kt

figure(5)
set(gcf, 'visible','off')

for i=0:numel(correlation_conductivity)-1
    subplot(numel(correlation_conductivity),3,3*i+1); contourf(T_check,P_check,squeeze(KT(:,:,1,i+1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
    xline(Tc,':','color','w')
    yline(Pc,':','color','w')
    hcb=colorbar;
    title(hcb,'kt');
    title('Conductivity 1-phase region:',' ')
    ylabel({correlation_conductivity{i+1};'Pressure [bar]'})
    xlabel('Temperature [K]')
    
    subplot(numel(correlation_conductivity),3,3*i+2); contourf(T_check,P_check,squeeze(KT(:,:,2,i+1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
    xline(Tc,':','color','w')
    yline(Pc,':','color','w')
    hcb=colorbar;
    title(hcb,'kt');
    title('Conductivity 2-phase region:','the biggest root')
    ylabel('Pressure [bar]')
    xlabel('Temperature [K]')
    
    subplot(numel(correlation_conductivity),3,3*i+3); contourf(T_check,P_check,squeeze(KT(:,:,3,i+1)), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; 
    xline(Tc,':','color','w')
    yline(Pc,':','color','w')
    hcb=colorbar;
    title(hcb,'kt');
    title('Conductivity 2-phase region:','the smallest root')
    ylabel('Pressure [bar]')
    xlabel('Temperature [K]')
end
print(gcf, '-dpdf', '-fillpage', 'KT.pdf'); close all

%}

%% Test RBF - one value

%{
Data = Z(:,:,1);

testpoint_P = 5;
testpoint_T = 5;

P_RBF = P_check(1:testpoint_P:end)/Pc;
T_RBF = T_check(1:testpoint_T:end)/Tc;

MU_RBF = [1e-3, 1e-2, 2e-2, 5e-2, 1e-1, 5e-1, 1e0];
APP = nan(numel(P_check), numel(T_check));

[X,Y] = meshgrid(P_RBF,T_RBF);
centers = [X(:), Y(:)];
distance = pdist2(centers,centers);

% Kernel Function
mu = MU_RBF(k);

GramMatrix = exp( - distance / mu );

% Weights
observations = Data(1:testpoint_P:end, 1:testpoint_T:end);
observations = reshape(observations.',[],1);

% Evaluate weights
%Weights = (GramMatrix'*GramMatrix)^(-1) * GramMatrix' * observations ;
Weights = GramMatrix \ observations ;

%% Evaluate the function at all the point

for i = 1:numel(P_check)
    for j = 1:numel(T_check)
        P = P_check(i)/Pc;
        T = T_check(j)/Tc;
        APP(i,j) = sum(Weights.*exp(- pdist2([P T], centers) / mu)');
    end
end
%}

%% Test RBF - loop over multiple values of mu

%{
Data = Z(:,:,1);

testpoint_P = 5;
testpoint_T = 5;

P_RBF = P_check(1:testpoint_P:end)/Pc;
T_RBF = T_check(1:testpoint_T:end)/Tc;

MU_RBF = [1e-3, 1e-2, 2e-2, 5e-2, 1e-1, 5e-1, 1e0];
APP = nan(numel(P_check), numel(T_check), numel(MU_RBF));

[X,Y] = meshgrid(P_RBF,T_RBF);
centers = [X(:), Y(:)];
distance = pdist2(centers,centers);

for k=1:numel(MU_RBF)

    if mod(k,round(numel(MU_RBF)/100)) == 0
        clc
        fprintf('%f %%\n',k/numel(MU_RBF)*100)
    end

    % Kernel Function
    mu = MU_RBF(k);

    GramMatrix = exp( - distance / mu );

    % Weights
    observations = Data(1:testpoint_P:end, 1:testpoint_T:end);
    observations = reshape(observations.',[],1);

    % Evaluate weights
    %Weights = (GramMatrix'*GramMatrix)^(-1) * GramMatrix' * observations ;
    Weights = GramMatrix \ observations ;

    %% Evaluate the function at all the point

    for i = 1:numel(P_check)

        for j = 1:numel(T_check)

            P = P_check(i)/Pc;
            T = T_check(j)/Tc;

            APP(i,j,k) = sum(Weights.*exp(- pdist2([P T], centers) / mu)');

        end

    end

    %max(max(abs(Data-APP)))

end

%%
for i=0:numel(MU_RBF)-1
    subplot(numel(MU_RBF),3,3*i+1);
    contourf(T_check,P_check,Data, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
    title('Original function')
    ylabel({['mu=',mat2str(MU_RBF(i+1))];'Pressure [bar]'})
    xlabel('Temperature [K]')

    subplot(numel(MU_RBF),3,3*i+2);
    hold on
    contourf(T_check,P_check,APP(:,:,i+1), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
    scatter(Y*Tc, X*Pc, 'k','SizeData',0.1)
    hold off
    title('RBF')
    ylabel('Pressure [bar]')
    xlabel('Temperature [K]')

    subplot(numel(MU_RBF),3,3*i+3);
    hold on
    contourf(T_check,P_check,Data-APP(:,:,i+1), 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
    scatter(Y*Tc, X*Pc, 'k','SizeData',0.1)
    hold off
    title('Difference')
    ylabel('Pressure [bar]')
    xlabel('Temperature [K]')
end
%}

%% Test RBF - optimze mu

mu = MX.sym('mu');

Data = RHO(:,:,1);

testpoint_P = 10;
testpoint_T = 10;

P_RBF = P_check(1:testpoint_P:end)/Pc;
T_RBF = T_check(1:testpoint_T:end)/Tc;

APP = MX(numel(P_check), numel(T_check));

[X,Y] = meshgrid(P_RBF,T_RBF);
centers = [X(:), Y(:)];
distance = pdist2(centers,centers);

% Kernel Function
GramMatrix = exp( - distance / mu );

% Weights
observations = Data(1:testpoint_P:end, 1:testpoint_T:end);
observations = reshape(observations.',[],1);

% Evaluate weights
disp(['Find the symbolic weights'])
%Weights = (GramMatrix'*GramMatrix)^(-1) * GramMatrix' * observations ;
tic
Weights = GramMatrix \ observations ;
toc

%% Evaluate the function at all the point

disp(['Create the matrix of data'])
for i = 1:numel(P_check)
    for j = 1:numel(T_check)
        P = P_check(i)/Pc;
        T = T_check(j)/Tc;
        APP(i,j) = sum(Weights.*exp(- pdist2([P T], centers) / mu)');
    end
end

%% Set the optimzation problem to find mu
disp(['Set the structure of the NLP'])
nlp = struct;            % NLP declaration
nlp.x = mu;              % decision vars

D = (Data-APP).^2;
MSE = sum(D(:))/numel(Data);
%MSE = norm(Data-APP,'fro')^2/numel(Data);
nlp.f = MSE;               % objective - mean squared error

disp(['Create the solver'])
% Create solver instance
F = nlpsol('F','ipopt',nlp);

disp(['Solve NLP'])
% Solve the problem using a guess
tic
res = F('x0',[1e-4]);
toc

%% Evaluate opt mu and corresponding weights
APP_opt = nan(numel(P_check), numel(T_check));

mu = full(res.x);
disp(['The optimal value of mu is ',num2str(mu)])

GramMatrix = exp( - distance / mu );
tic
Weights = GramMatrix \ observations ;
toc

%% Approximate the dataset
for i = 1:numel(P_check)
    for j = 1:numel(T_check)
        P = P_check(i)/Pc;
        T = T_check(j)/Tc;
        APP_opt(i,j) = sum(Weights.*exp(- pdist2([P T], centers) / mu)');
    end
end

%% Plot the results with opt mu and compare to the original data
subplot(1,3,1);
contourf(T_check,P_check, Data, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
title('Original function')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,2);
hold on
contourf(T_check,P_check,APP_opt, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
scatter(Y*Tc, X*Pc, 'k')
hold off
title('RBF')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')

subplot(1,3,3);
hold on
contourf(T_check,P_check,Data-APP_opt, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar
scatter(Y*Tc, X*Pc, 'k')
hold off
title('Difference')
ylabel('Pressure [bar]')
xlabel('Temperature [K]')