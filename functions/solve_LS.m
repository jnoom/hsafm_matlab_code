function [Fhhat,xhat] = solve_LS(Liy,LiF,n,Fhhat1,err)
NT        = length(Fhhat1(1,:));
[row,col] = find(max(0,abs(Fhhat1)-err));
xF        = zeros(length(LiF(1,:)),NT);

K                                = true(size(Fhhat1));
K(sub2ind(size(Fhhat1),row,col)) = false;

for i = 1:NT
    K2       = logical(1 - [kron(K(:,i),[false(n,1);1]);false(n,1)]);  % Determine nonzero components
    xF(K2,i) = LiF(:,K2) \ Liy(:,i);                                   % solve the (reduced) LS problem
end

Fhhat = xF(n+1:n+1:end,:);                % Extract estimated t/s-interactions from xF
for i = 1:n
    xhat(i:n:n*NT,:) = xF(i:n+1:end,:)';  % Extract estimated system states from xF
end
end