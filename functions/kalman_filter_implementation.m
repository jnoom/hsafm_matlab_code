function theta = kalman_filter_implementation(A,C,u,y,Rgamma,x0,P0,Q,kT,k1) 
n      = length(C(1,:));
P      = P0;
x      = [x0,zeros(n,length(y)-1)];
ye     = zeros(length(y),1);
errcov = zeros(length(y),1); 

for i = kT(n+1):length(y) 
  % Measurement update
  Mn     = P*C(i,:)'/(C(i,:)*P*C(i,:)'+Rgamma);
  x(:,i) = x(:,i) + Mn*(y(i)-C(i,:)*x(:,i));    % x[n|n]
  P      = (eye(n)-Mn*C(i,:))*P;                % P[n|n]

  ye(i)        = C(i,:)*x(:,i);
  errcov(i)    = C(i,:)*P*C(i,:)';
  theta(1:n,i) = x(:,i);
  
  % Time update
  if i-kT(k1(i))<1
    x(:,i+1) = A(:,:,2)*x(:,i) + u(i);             % x[n+1|n]
    P        = A(:,:,2)*P*A(:,:,2)' + Q(:,:,2);    % P[n+1|n] 
  else
    x(:,i+1) = A(:,:,1)*x(:,i) + u(i);             % x[n+1|n]
    P        = A(:,:,1)*P*A(:,:,1)' + Q(:,:,1);    % P[n+1|n]
  end
end