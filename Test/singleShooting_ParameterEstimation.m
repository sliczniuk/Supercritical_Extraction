function [xout, kout] = singleShooting_ParameterEstimation(OCP, x0, u, yd, k0)

%addpath('C:\Users\slicz\Desktop\Aalto\Matlab\SFE\casadi-windows-matlabR2016a-v3.5.1');
%addpath('\\home.org.aalto.fi\sliczno1\data\Documents\casadi-windows-matlabR2016a-v3.5.1');
import casadi.*

OCP_solver = casadi.Opti();
%OCP_solver.solver('ipopt')

% http://www.diva-portal.se/smash/get/diva2:956377/FULLTEXT01.pdf
nlp_opts = struct;
nlp_opts.ipopt.max_iter = 100;
%nlp_opts.ipopt.acceptable_iter = 50;
%nlp_opts.ipopt.acceptable_tol = 1e-6;
%nlp_opts.ipopt.tol = 1e-7;
ocp_opts = {'nlp_opts', nlp_opts};

OCP_solver.solver('ipopt',nlp_opts)

K = OCP_solver.variable(OCP.Nk);
X = [MX(x0) zeros(OCP.Nx,OCP.N)];
J = 0;
for j=1:OCP.N
    %J=J+OCP.L(U(:,j));
    X(:,j+1)=OCP.F(X(:,j),u(:,j),K);
end
J = J + OCP.Lf(X,yd);

% for nx=1:OCP.Nx
%     OCP_solver.subject_to(OCP.x_lu(nx,1)<=X(nx,:)<= OCP.x_lu(nx,2));
% end
% 
% for nk=1:OCP.Nk
%     OCP_solver.subject_to(OCP.k_lu(nk,1)<=K(nk,:)<= OCP.k_lu(nk,2));
% end

OCP_solver.minimize(J);
OCP_solver.set_initial(K, k0);

try
    sol = OCP_solver.solve();
    kout = sol.value(K);
catch
    kout = OCP_solver.debug.value(K);
end

xout = 1;
end