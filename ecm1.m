function [x,y,phi] = ecm1(xp,u,up,par,ocv)
    N = length(ocv)-1;
    Ts = 1;
    R0 = par(1);
    R1 = par(2);
    C1 = par(3);
    Qn = 3.08;

    A = [1-Ts/(R1*C1) 0;0 1];
    Bu = [Ts/C1;-Ts/(3600*Qn)];
    C =[-1 0];
    Du = [-R0];

    x = A*xp + Bu*up;
    phi = ocv*(x(2,1).^[N:-1:0])';
    y = C*x+Du*u+phi;
end