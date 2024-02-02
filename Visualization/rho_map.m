startup;
delete(gcp('nocreate'));
%p = Pushbullet(pushbullet_api);

%addpath('C:\Dev\casadi-3.6.3-windows64-matlab2018b');
addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-3.6.3-windows64-matlab2018b');
import casadi.*

%%

Parameters_table        = readtable('Parameters.csv') ;             % Table with prameters
Parameters              = num2cell(Parameters_table{:,3});          % Parameters within the model + (m_max), m_ratio, sigma

T                       = 300:1:330 ;
P                       = 200:5:300 ;
[X,Y]                   = meshgrid(T,P);
data                    = [X(:), Y(:)];
RHO                     = [];

for i = 1:length(data)
    Z                   = Compressibility( data(i,1), data(i,2),         Parameters );
    rho                 = full(rhoPB_Comp( data(i,1), data(i,2), Z,      Parameters ));
    RHO                 = [RHO, rho];
end

%%
excel_file              = 'output_T45.xls';
data_opt                = readcell(excel_file,'sheet',33);
TT_out                  = [cell2mat(data_opt(2:151,3))];
P_opt                   = [cell2mat(data_opt(2:151,4))];

%%
s                       = pcolor(X,Y,reshape(RHO,length(P),[]));
s.EdgeColor             = 'none';

hold on
plot(TT_out,P_opt,'ok')
hold off