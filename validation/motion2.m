%% define Eulerian motion function
function [v, dv] = motion2(xIJ, t)
    
    gamma1 = 0.6;
    gamma2 = 1. - gamma1;

    Omega1 = [0., -1.; 1., 0.]*pi/4.;
    Omega2 = [0., -1.; 1., 0.]*pi/40.;

    X1 = [0.5, 0.5]';
    X2 = [0.8, 0.8]';

    x0 = gamma1 * X1 + gamma2 * X2;

    Q1 = expm(t*Omega1);
    Q2 = expm(t*Omega2);

    R = gamma1 * Q1 + gamma2 * Q2;
    Rinv = inv(R);

    S1 = Q1 * Rinv;
    S2 = Q2 * Rinv;

    x1 = Q1 * X1;
    x2 = Q2 * X2;

    xTilde = gamma1 * x1 + gamma2 * x2 - x0;

    v1 = Omega1 * (S1 * (xIJ + xTilde) - x1);
    v2 = Omega2 * (S2 * (xIJ + xTilde) - x2);
    v = gamma1 * v1 + gamma2 * v2;

    gradv1 = Omega1 * S1;
    gradv2 = Omega2 * S2;
    gradv  = gamma1 * gradv1 + gamma2 * gradv2;

    dv = gamma1 * (Omega1 * v1) + gamma2 * (Omega2 * v2) - gradv * v;

end
