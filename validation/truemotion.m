function x = truemotion(X, t)
    
    omega = [ 0, -pi; pi, 0 ];
    V0 = [ 0.1; 0.1 ];
    X0 = [ 0.5; 0.5 ];

    Q = expm(t * omega);

    x  = Q * (X - X0) + X0 + V0*t;

end
