%% define Eulerian motion function
function [v, dv] = motion(x, t)
    
    omega = [ 0, -pi; pi, 0 ];
    V0 = [ 0.1; 0.1 ];
    X0 = [ 0.5; 0.5 ];

    v  = omega * (x - X0 - V0 * t) + V0;
    dv = -omega * V0;
end
