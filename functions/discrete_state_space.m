function [xnew,y] = discrete_state_space(A,B,C,D,xold,u)
xnew = A*xold + B*u;
y    = C*xold + D*u;