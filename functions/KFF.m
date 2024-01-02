function Fhhat_KFF = KFF(yc, uc, sys, Td, Rdelta, Rgamma,n,delay)
fprintf('Applying KFF... ')
A=sys.A; B=sys.B; C=sys.C; D=sys.D; ts=sys.Ts;
l = length(yc);

%% Find time instances for tip-sample interaction
j=1;
for i=1:l-1
    if uc(i) < 0 && uc(i+1) >= 0
        k(j) = i;
        j    = j+1;
    end
end
kT = k-delay;

delta = 0;
j     = n+1;
for i=kT(n):l
    if i ~= kT(j)
        k1(i) = find(kT==i-delta);
        delta = delta+1;
    else
        delta = 0;
        j     = min(j+1,length(kT));
        k1(i) = find(kT==i-delta);
        delta = delta+1;
    end
end

%% Construct Al and Cl
fprintf('constructing matrix Cl')
for i=1:kT(2*n)+200
    CAB(i) = C*A^(i-1)*B;
end
for i=kT(n):l
    if i-kT(k1(i))<1
        Cl(i,1) = 0;
        for k=2:n
            Cl(i,k) = CAB(i-kT(k1(i)-k+1));
        end
    else
        for k=1:n
            Cl(i,k) = CAB(i-kT(k1(i)-k+1));
        end
    end
end

%% Kalman filter implementation
fprintf(repmat('\b', 1, 22)); fprintf('calculating Fhhat')
Al(:,:,1) = eye(n);
Al(:,:,2) = [1 zeros(1,n-2), 0;
              eye(n-1),      zeros(n-1,1)];
Q(:,:,1)  = zeros(n);
Q(:,:,2)  = [Rdelta, zeros(1,n-1);
             zeros(n-1,n)];
u  = zeros(l,1);
x0 = zeros(length(Cl(1,:)),1);
P0 = eye(length(Cl(1,:)));

for k=kT(find(kT>2*Td/ts,1)):l
    fact(k) = 0;
    for i=k-2*round(Td/ts):k-1
        fact(k) = fact(k) + CAB(k-i)*uc(i);
    end
    Psi(k) = yc(k) - fact(k);
end

theta = kalman_filter_implementation(Al,Cl,u,Psi,Rgamma,x0,P0,Q,kT,k1);

Fhhat_KFF = zeros(l,1);
for k=kT(n):l
    if k > kT(end-n)
        if k > kT(end-n+1)
            if n > 2
                if k > kT(end-n+2)
                    if n > 3
                        if k > kT(end-n+3)
                            if n > 4
                                if k > kT(end-n+4)
                                    Fhhat_KFF(k) = 0;
                                else
                                    Fhhat_KFF(kT(k1(k))) = theta(n-4,kT(k1(k)+n-4)-1);
                                end
                            else
                                Fhhat_KFF(k) = 0;
                            end
                        else
                            Fhhat_KFF(kT(k1(k))) = theta(n-3,kT(k1(k)+n-3)-1);
                        end
                    else
                        Fhhat_KFF(k) = 0;
                    end
                else
                    Fhhat_KFF(kT(k1(k))) = theta(n-2,kT(k1(k)+n-2)-1);
                end
            else
                Fhhat_KFF(k) = 0;
            end
        else
            Fhhat_KFF(kT(k1(k))) = theta(n-1,kT(k1(k)+n-1)-1);
        end
    else
    Fhhat_KFF(kT(k1(k))) = theta(n,kT(k1(k)+n)-1);
    end
end
fprintf(repmat('\b', 1, 17));disp('Done!')
end