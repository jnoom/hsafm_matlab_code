function [LiF,Liy] = fista_matrices(u,y,sys,Rnu,Reta,n,T)
A = sys.A; B = sys.B; C = sys.C; G = sys.B;

Rni = inv(sqrtm(Rnu));
RnC = Rni*C;
Rei = eye(n)*inv(sqrtm(Reta));
ReA = Rei*A;
ReG = Rei*G;

ny = length(C(:,1));
nf = length(G(1,:));

LiF1 = [RnC, zeros(ny,nf), zeros(ny,n);
        -ReA, -ReG, Rei];
Liy  = [Rni*y(1,:);Rei*B*u(1,:)];
LiF  = LiF1;
for l = 2:T
    LiF((ny+n)*(l-1)+1:(ny+n)*l, (n+nf)*(l-1)+1:(n+nf)*l+n) = LiF1;
    Liy((ny+n)*(l-1)+1:(ny+n)*l,:) = [Rni*y(l,:);Rei*B*u(l,:)];
end
LiF = LiF(1:end-n,1:end-n-nf);
Liy = Liy(1:end-n,:);

end

