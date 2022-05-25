function xdot = test_model(x, u, k, theta)
    
    Ca = x(1);
    Cb = x(2);
    
    k1 = k(1);
    k2 = k(2);

    xdot = [-k1*Ca;
            k1*Ca - k2*Cb;
            k2*Cb + u];

end