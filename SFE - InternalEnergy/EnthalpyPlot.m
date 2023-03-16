clc, close all
clear all

%addpath('C:\dev\casadi-windows-matlabR2016a-v3.5.2');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

Parameters_table        = readtable('Parameters.csv') ;        % Table with prameters
Parameters_cell         = table2cell(Parameters_table(:,3));

%%
T_set               = linspace(34,100)+273;
P_set               = linspace(74,300,200);

data                = [combvec(T_set,P_set)'];

h_data              = nan(numel(T_set)*numel(P_set),1);

parfor i=1:length(data)
    T               = data(i,1);
    P               = data(i,2);

    Z               = Compressibility( T, P,         Parameters_cell );
    rho             = rhoPB_Comp(      T, P, Z,      Parameters_cell );
    h               = SpecificEnthalpy(T, P, Z, rho, Parameters_cell );
    h_data(i)       = h;        
end

%%
[X,Y] = meshgrid(T_set,P_set);
Z     = reshape(h_data,numel(T_set),numel(P_set))'*1e-3;

%%
figure(1);
set(gcf,'PaperOrientation','landscape', 'visible','on')

contourf(X,Y,Z, 50, 'EdgeColor', 'none'); colormap jet; axis square tight; colorbar;

%xline((Parameters_cell{10}),':','color','w')
%yline(Parameters_cell{11},':','color','w')

hcb=colorbar;
title(hcb,'Enthalpy~[kJ/kg]','Interpreter','latex');

ylabel('Pressure [bar]')
xlabel('Temperature [K]')

print(gcf, '-dpdf', '-fillpage', 'Enthalpy.pdf'); close all