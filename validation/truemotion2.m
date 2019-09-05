function x = truemotion2(X, t)
    
    gamma1 = 0.6;
    gamma2 = 1. - gamma1;

    Omega1 = [0., -1.; 1., 0.]*pi/4.;
    Omega2 = [0., -1.; 1., 0.]*pi/40.;

    X1 = [0.5, 0.5]';
    X2 = [0.8, 0.8]';

    Q1 = expm(t*Omega1);
    Q2 = expm(t*Omega2);

    x1 = Q1 * (X - X1) + X1;
    x2 = Q2 * (X - X2) + X2;
    x = gamma1 * x1 + gamma2 * x2;

end
