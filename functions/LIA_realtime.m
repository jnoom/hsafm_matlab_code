function [A_LIA, phi_LIA, xnewx, xnewy] = LIA_realtime(yc,ts,f_drive,g,Ab,Bb,Cb,Db,xoldx,xoldy)
l   = length(yc);

Vi  = g*yc(l);
vrx = cos(2*pi*f_drive*(l-1)*ts)';
vry = sin(2*pi*f_drive*(l-1)*ts)';

Vmx = Vi*vrx;
Vmy = Vi*vry;

[xnewx,Voutx] = discrete_state_space(Ab,Bb,Cb,Db,xoldx,Vmx);
[xnewy,Vouty] = discrete_state_space(Ab,Bb,Cb,Db,xoldy,Vmy);

A_LIA   = (2/g) * sqrt(Voutx.^2 + Vouty.^2);
phi_LIA = atan(Vouty./Voutx);
end