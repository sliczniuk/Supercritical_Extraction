function k = HeatConductivity_Comp(T,P,Z,rho, theta, correlation)

    switch correlation
        case 'Amooey'
            A1 = theta{28}; % -105.161;
            A2 = theta{29}; % 0.9007;
            A3 = theta{30}; % 0.0007;
            A4 = theta{31}; % 3.50E-15;
            A5 = theta{32}; % 3.76E-10;
            A6 = theta{33}; % 0.7500;
            A7 = theta{34}; % 0.0017;
    
            k = ( A1 + A2.*rho + A3.*(rho.^2) + A4.*(rho.^3).*(T.^3) + A5.*(rho.^4) + A6.*T + A7.*(T.^2) ) ./ (T.^(1/2));
    
        case 'Bahadori'
            P = P/10;
            a =  2.51177    - (4.61299e3)/T + (1.5604e6)  /(T^2) - (1.64868e8 )/(T^3);
            b = -6.78436e2  + (5.94729e5)/T - (1.81369e8) /(T^2) + (1.86064e10)/(T^3);
            c =  2.064898e4 - (1.99667e7)/T + (6.42367e9) /(T^2) - (6.8022e11 )/(T^3);
            d = -1.09504e5  + (1.08783e8)/T - (3.57549e10)/(T^2) + (3.855e12  )/(T^3);
    
            Lnk = a + b/P + c/(P^2) + d/(P^3);
    
            k = exp(Lnk)*10^3;
    
        case 'Jarrahian'
            P = P/10;
            A1= 1.49288e+1;
            A2= 2.62541e-3;
            A3= 8.77805e-6;
            A4=-5.11425;
            A5= 4.37711e-1;
            A6= 2.11405e-5;
            A7=-4.73036e-1;
            A8= 7.36636e-2;
            A9=-3.76340e-3;
    
            k = ( A1 + A2*P + A3*(P^2) + A4*log(T) + A5*(log(T)^2) ) / ( 1 + A6*P + A7*log(T) + A8*(log(T)^2) + A9*(log(T)^3) );
    
        case 'Rostami'
            P = P/10;
            if P < 20
                k = 0.0936*T - 0.448*P + 0.0739*rho - 0.244*log(rho) - 10.8*P/T + 0.00753*P^2 + 1.85e-5*rho^2 + 94.7*rho/(T*P) - 16.7;
            else
                k = 0.0575*T + 0.0151*P + 0.0372*rho + 4.69e-5*log(log(rho)) - 0.00695*log(P) + 1.41e-5*T^2 + 7.2e-8*rho^3 + 1.78;
            end
    
        case 'Rostamian'
            A1  =  -29.9717451505165e+015;
            A2  =    9.65637447009372e+018;
            A3  =  -13.8288944829492e+021;
            A4  =  -21.1152877719961e+012;
            A5  =    9.26006733304733e+000;
            A6  =  -30.7171646680127e+006;
            A7  = -408.256276723566e+015;
            A8  =  130.491020289031e+018;
            A9  =  -13.4237924607890e+021;
            A10 =   60.9547298940653e+009;
    
            k = ( A1 + A2/T + A3/(T^2.5) + A4*rho + A5*(rho^2.5)*(T^2.5) + A6*(rho^3) ) / ( 1 + A7/T + A8/(T^2.5) + A9/(T^3) + A10*sqrt(rho*T) );

        case 'Huber'
            Tc =  304.1282;
            rho_c = 467.6;

            Tr = T/Tc;

            Lk = [1.51874307e-2, 2.80674040e-2, 2.28564190e-2, -7.41624210e-3];

            B1 = [1.00128e-2, 5.60488e-2, -8.11620e-2, 6.24337e-2, -2.06336e-2,  2.53248e-3];
            B2 = [4.30829e-3, -3.58563e-2, 6.71480e-2, -5.22855e-2, 1.74571e-2,  -1.96414e-3];

            Delta_Tc = Tr - 1;
            Delta_rho_c = rho/rho_c - 1;

            k0 = sqrt(Tr) / ( Lk(1)/(Tr^0) + Lk(2)/(Tr^1) + Lk(3)/(Tr^2) + Lk(4)/(Tr^3) );

            Delta_k = 0;
            for i=1:6
                Delta_k = (B1(i) + B2(i)*Tr) * (rho/rho_c)^i;
            end

            Delta_k_c = (-17.47-44.99 * Delta_Tc) / ( 0.8563 - exp(8.865*Delta_Tc + 4.16*Delta_rho_c^2 + 2.302*Delta_Tc*Delta_rho_c - Delta_rho_c^3) - 0.4503*Delta_rho_c - 7.197*Delta_Tc );

            k = k0 + Delta_k + Delta_k_c;

        otherwise
            disp('No such a correlation')
            
    end

end