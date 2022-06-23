%% Newtonian Motion KF
close all
clear all
clc

A = [1 .1; 0 1];
H = [1 0; 0 1];
G = [0 1];
B = [.01/2; .1];
u = [9.8];
R = [2/.1];
Q = ([1 0;0 1]*.1) + (((.1)^2)*[0 .5; .5 0]) + (((.1)^3)*[1/3 0; 0 0]);
x = [0; 0];
P = [1 0; 0 1];
z = [0;0];
for i=1:100
    r = rand
    z(1,1) = z(1,1) + .44 + (.1*r)^2;
    z(2,1) = z(2,1) + .98;
    Pp = A*P*A' + G*Q*G';
    P = [Pp^-1 + H'*(R^-1)*H]^-1;
    xp = A*x + B*u
    z
    zr = z - H*xp
    x = xp + P*H'*(R^-1)*(zr)
end
Pp
P
x
