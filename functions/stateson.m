% Inputs:
% -u2           Cantilever input
% -y2           Measured cantilever deflection
% -sys          Identified (state-space) system
% -n            Order of identified system
% -Rnu          Estimated variance of measurement noise
% -Reta         Estimated variance of process noise
% -Q            Estimated variance of pulses in tip-sample interaction
% -T            Number of measurements per STATESON evaluation
% -lambdafactor Factor for calculating lambda
% -epsfactor    Factor for epsilon
% -err          Minimum of estimated pulses
%
% Outputs:
% -Fhhat1       Estimated tip-sample interaction
% -yhat1        Estimated cantilever deflection

function [Fhhat1,yhat1] =...
    stateson(u2,y2,sys,n,Rnu,Reta,Q,T,lambdafactor,epsfactor,err)
i=1;
while i*T < length(u2)
    u(:,i) = u2((i-1)*T+1:i*T);
    y(:,i) = y2((i-1)*T+1:i*T);
    i      = i+1;
end
NT = i;       % Number of data segments
fprintf('Performing STATESON in %d steps... ',NT)

u(:,i) = [u2((i-1)*T+1:end); zeros(T-length(u2((i-1)*T+1:end)),1)];
y(:,i) = [y2((i-1)*T+1:end); zeros(T-length(y2((i-1)*T+1:end)),1)];

[LiF,Liy] = fista_matrices(u,y,sys,Rnu,Reta,n,T);  % Construct matrices for FISTA evaluation
  
alpha = ones(T-1,1);
eps   = epsfactor*ones(T-1,1);

Fhhat = solve_lasso(Liy,LiF,Q,alpha,lambdafactor,n);  % lasso problem
for i=1:NT
    alpha(:,i) = 1./(eps+norms(inv(sqrtm(Q))*Fhhat(:,i),1,2));
end

[Fhhat,X] = solve_LS(Liy,LiF,n,Fhhat,err);  % (reduced) least-squares problem

for i=1:NT
    yhat(:,i) = sys.C*X(n*(i-1)+1:n*i,:);
end

disp(' ')
Fhhat1 = vec([Fhhat;zeros(1,NT)]);
yhat1 = vec(yhat(1:T,:));
end