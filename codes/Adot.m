function Adot = Adot(A, t)

Adot = [1 sin(A(1))*tan(A(2)) cos(A(1))*tan(A(2)); 0 cos(A(1)) -sin(A(1)); 0 sin(A(1))/cos(A(2)) cos(A(1))/cos(A(2))]*...
    [sin(0.1*t); 0.01; cos(0.1*t) ]*20*pi/180