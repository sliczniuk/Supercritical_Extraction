startup;

import casadi.*

Parameters_table        = readtable('Parameters.csv') ;        % Table with prameters
Parameters_cell         = table2cell(Parameters_table(:,3));

%%
T_set               = linspace(34, 137, 100)+273;
P_set               = linspace(74, 300, 100);

[p,q]               = meshgrid(T_set,P_set);
data                = [p(:) q(:)];

h_data              = nan(numel(T_set)*numel(P_set),1);
Z_data              = nan(numel(T_set)*numel(P_set),1);
rho_data            = nan(numel(T_set)*numel(P_set),1);
cp_data             = nan(numel(T_set)*numel(P_set),1);
mu_data             = nan(numel(T_set)*numel(P_set),1);
k_data              = nan(numel(T_set)*numel(P_set),1);

numIterations       = length(data);

f                   = waitbar(0, 'Starting');

for i=1:numIterations
    T               = data(i,1);
    P               = data(i,2);

    Z               = Compressibility( T, P,         Parameters_cell );
    rho             = rhoPB_Comp(      T, P, Z,      Parameters_cell );
    h               = SpecificEnthalpy(T, P, Z, rho, Parameters_cell );
    Cp              = SpecificHeatComp(T, P, Z, rho, Parameters_cell );
    mu              = Viscosity(T,rho);
    k               = HeatConductivity_Comp(T, rho);

    Z_data(i)       = Z;
    rho_data(i)     = rho;
    h_data(i)       = h;     
    cp_data(i)      = Cp;     
    mu_data(i)      = mu;     
    k_data(i)       = k;     

    waitbar(i/numIterations, f, sprintf('Progress: %d', floor(i/numIterations*100)));
    pause(0.1);
end

close(f)

%%
T_set = T_set - 273;
[X,Y] = meshgrid(T_set,P_set);

%%
Z     = reshape(h_data,numel(T_set),numel(P_set));

figure(1);

set(gcf,'PaperOrientation','landscape', 'visible','on')

[C,h] = contourf(X,Y,Z, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar;

%clabel(C,h)
%title('$\textbf{Enthalpy,~1-phase~region}$')

%xline((Parameters_cell{10}),':','color','w')
%yline(Parameters_cell{11},':','color','w')

hcb=colorbar;
title(hcb,'$h~[kJ/kg]$','Interpreter','latex');

ylabel('Pressure~[bar]')
xlabel('Temperature~[C]')

fontsize(36,"points")
print(gcf, '-dpdf', '-fillpage', 'Enthalpy.pdf'); close all

%%
Z     = reshape(Z_data,numel(T_set),numel(P_set));

figure(1);

set(gcf,'PaperOrientation','landscape', 'visible','on')
set(gca,'FontName','Segoe UI')

[C,h] = contourf(X,Y,Z, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar;

%clabel(C,h)
%title('$\textbf{Compressibility,~1-phase~region}$')

%xline((Parameters_cell{10}),':','color','w')
%yline(Parameters_cell{11},':','color','w')

hcb=colorbar;
title(hcb,'$Z~[-]$','Interpreter','latex');

ylabel('Pressure~[bar]')
xlabel('Temperature~[C]')

fontsize(36,"points")
print(gcf, '-dpdf', '-fillpage', 'Compressibility.pdf'); close all

%%
Z     = reshape(rho_data,numel(T_set),numel(P_set));

figure(1);

set(gcf,'PaperOrientation','landscape', 'visible','on')

[C,h] = contourf(X,Y,Z, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar;

%clabel(C,h)
%title('$\textbf{Density,~1-phase~region}$')

%xline((Parameters_cell{10}),':','color','w')
%yline(Parameters_cell{11},':','color','w')

hcb=colorbar;
title(hcb,'$\rho~[kg/m^3]$','Interpreter','latex');

ylabel('Pressure~[bar]')
xlabel('Temperature~[C]')

fontsize(36,"points")
print(gcf, '-dpdf', '-fillpage', 'RHO.pdf'); close all

%%
Z     = log10(reshape(cp_data,numel(T_set),numel(P_set)));

figure(1);

set(gcf,'PaperOrientation','landscape', 'visible','on')

[C,h] = contourf(X,Y,Z, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar;

%clabel(C,h)
%title('$\textbf{Specific~heat,~1-phase~region}$')

%xline((Parameters_cell{10}),':','color','w')
%yline(Parameters_cell{11},':','color','w')

hcb=colorbar;
title(hcb,'$log_{10}C_p~[kJ/kg/K]$','Interpreter','latex');

ylabel('Pressure~[bar]')
xlabel('Temperature~[C]')

fontsize(36,"points")
print(gcf, '-dpdf', '-fillpage', 'CP.pdf'); close all

%%
Z     = reshape(mu_data*1e5,numel(T_set),numel(P_set));

figure(1);

set(gcf,'PaperOrientation','landscape', 'visible','on')

[C,h] = contourf(X,Y,Z, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar;

%clabel(C,h)
%title('$\textbf{Viscosity,~1-phase~region}$')

%xline((Parameters_cell{10}),':','color','w')
%yline(Parameters_cell{11},':','color','w')

hcb=colorbar;
title(hcb,'$\mu~[Pa~s] \times 10^{-5}$','Interpreter','latex');

ylabel('Pressure~[bar]')
xlabel('Temperature~[C]')

fontsize(36,"points")
print(gcf, '-dpdf', '-fillpage', 'MU.pdf'); close all

%%
Z     = reshape(k_data,numel(T_set),numel(P_set));

figure(1);

set(gcf,'PaperOrientation','landscape', 'visible','on')

[C,h] = contourf(X,Y,Z, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar;

%clabel(C,h)
%title('$\textbf{Heat conductivity,~1-phase~region}$')

%xline((Parameters_cell{10}),':','color','w')
%yline(Parameters_cell{11},':','color','w')

hcb=colorbar;
title(hcb,'$k~[(W/mK)]\times 10^{-3}$','Interpreter','latex');

ylabel('Pressure~[bar]')
xlabel('Temperature~[C]')

fontsize(36,"points")
print(gcf, '-dpdf', '-fillpage', 'KT.pdf'); close all