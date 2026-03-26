function [x,y,phi] = ecm2(xp,u,up,par,ocv)
    N = length(ocv)-1;
    Ts = 1;
    R0 = par(1);
    R1 = par(2);
    R2 = par(3);
    C1 = par(4);
    C2 = par(5);
    Qn = 3.08;

    A = [1-Ts/(R1*C1) 0 0;0 1-Ts/(R2*C2) 0;0 0 1];
    Bu = [Ts/C1;Ts/C2;-Ts/(3600*Qn)];
    C =[-1 -1 0];
    Du = [-R0];

    x = A*xp + Bu*up;
    phi = ocv*(x(3,1).^[N:-1:0])';
    y = C*x+Du*u+phi;
end