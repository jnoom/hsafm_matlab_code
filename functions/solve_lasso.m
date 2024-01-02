function Fhhat = solve_lasso(Liy,LiF,Q,w,lambda,n)
opts.pos          = false;
opts.lambda       = 1/(sqrtm(Q)).*[kron(w,[zeros(n,1);1]);zeros(n,length(w(1,:)))]*lambda;
opts.backtracking = false;

xF    = fista_lasso(Liy, LiF, [], opts);  % from FISTA package
Fhhat = xF(n+1:n+1:end,:);
end