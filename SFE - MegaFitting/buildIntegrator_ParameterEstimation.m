function F = buildIntegrator_ParameterEstimation(varargin)
    
    import casadi.*

    f = varargin{1};
    d = varargin{2};
    
    % Verify if additional arguments were provided
    
    options = struct;
    
    if(numel(varargin) > 2)
        options.tf = varargin{3};
    else
        options.tf = 1;     
    end
    
    if(numel(varargin) > 3)
        method = varargin{4};     
    else                    
        method = 'cvodes';     
    end
    
    x = MX.sym('x',d(1));
    u = MX.sym('u',d(2));
    k = MX.sym('k',d(3));
    p = [u];
    
    DAE      = struct('x', x, 'p', [u; k], 'ode', f(x,u, k));
    args     = { x , u , k };
    args_str = {'x','u','k'};
    
    %options.nonlinear_solver_iteration = 'functional';
    %options.newton_scheme = 'tfqmr';
    
    F     = integrator('F', method, DAE, options);
    F_res = F('x0', x, 'p', [u; k]);
    
    F     = Function('F', args, {F_res.xf}, args_str, {'x_next'});
    
end