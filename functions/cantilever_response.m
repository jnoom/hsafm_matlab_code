% Inputs: see "hsafm_stateson.m"
%
% Outputs:
% -yc2          Measured cantilever deflection
% -xpos2        Image coordinates
% -ypos2        Image coordinates
% -z            Scanner height
% -hhat_miny    Estimated sample height using MRM
% -hhat_uz      Estimated sample height using conventional method
function [yc2, xpos2, ypos2, z, hhat_miny, hhat_uz] = ...
    cantilever_response(sys2,x0,uc2,Reta,Rnu,h,c,scansize,scanspeed,lines,f_drive,fcut,controller,zpiezo,Ampref)
%% Initialization
l2           = length(uc2);  
samplefactor = 50;                      % Upsampling factor for simulating cantilever dynamics (simulation sampling time << measurement sampling time)
ts2          = sys2.Ts;                 % Measurement sampling time

CZ=controller*zpiezo; Acz=CZ.A; Bcz=CZ.B; Ccz=CZ.C; Dcz=CZ.D;        % Controller & z-piezo dynamics
Ac=controller.A; Bc=controller.B; Cc=controller.C; Dc=controller.D;  % Controller dynamics

nbut=4; [Ab,Bb,Cb,Db] = butter(nbut,fcut*ts2*2,'low');  % Filter for lock-in amplifier
xnewx=zeros(size(Bb)); xnewy=xnewx;     % Initial condition of lock-in amplifier

sys  = c2d(d2c(sys2),ts2/samplefactor); % Discretized system for simulation of cantilever dynamics
A=sys.A; B=sys.B; C=sys.C; D=sys.D; ts=sys.Ts;  
eta2 = sqrt(Reta) * randn(l2,1);        % Measurement noise
nu2  = sqrt(Rnu) * randn(l2,1);         % Process noise
ts3  = 1/f_drive;                       % Cantilever oscillation period

for i=1:round(length(uc2)*ts2/ts)       % Convert to higher sampling rate
    uc(i,1)  = uc2(ceil(i*ts/ts2));
    eta(i,1) = eta2(ceil(i*ts/ts2));
end
l   = length(uc);
n   = length(B);
ncz = length(Bcz);
nc  = length(Bc);

in    = scansize*2/scanspeed/ts;  % Determine image coordinates
in2   = scansize*2/scanspeed/ts2;
xpos  = max(1,ceil((sawtooth([1:l]/in*2*pi,0.5)/2+0.5) * length(h(1,:))));
ypos  = floor(linspace(1,lines+1,l));
xpos2 = max(1,ceil((sawtooth([1:l2]/in2*2*pi,0.5)/2+0.5) * length(h(1,:))));
ypos2 = floor(linspace(1,lines+1,l2));

x            = zeros(n,l+1);  % Initial conditions
x(1:n,1)     = x0;
xc           = zeros(nc,l2+1);
xc(1:nc,1)   = zeros(nc,1);
xcz          = zeros(ncz,l2+1);
xcz(1:ncz,1) = zeros(ncz,1);
yc           = zeros(l,1); 
F            = zeros(l,1);
prog         = 0;

j=2; k=1; k1=1; k2=round(ts3/ts2)+1; Amp2(1)=0; hz=h;

fprintf('Simulating cantilever...  0%%');

for i = 1:l    
    %% Cantilever deflection
    x(1:n,i+1) = A*x(1:n,i) + B*(uc(i)+eta(i));
    yc(i)      = C*x(1:n,i) + D*uc(i);          % Cantilever deflection in case of no impact
    
    if yc(i) < hz(min(lines,ypos(i)),xpos(i))   % Calculate t/s-interaction (piece-wise linear)
        F(i) = c * (hz(min(lines,ypos(i)),xpos(i)) - yc(i));
    else
        F(i) = 0;
    end
    
    x(1:n,i+1) = A*x(1:n,i) + B*(uc(i)+eta(i)+F(i));
    yc(i)      = C*x(1:n,i) + D*uc(i);               % Real cantilever deflection

    %% Lock-in amplifier & control of scanner height
    if i == k1
        yc2(k,1) = yc(i)+nu2(k);    % Measured cantilever deflection
        if k > 1
            [Amp1(k),~,xnewx,xnewy] = LIA_realtime(yc2(1:k),ts2,f_drive,7,Ab,Bb,Cb,Db,xnewx,xnewy);
            if k == k2              % One LIA measurement per oscillation period
                Amp2(k) = Amp1(k);
                k2      = k2+round(ts3/ts2);
            else
                Amp2(k) = Amp2(k-1);
            end
            xc(1:nc,k+1)   = Ac*xc(1:nc,k) + Bc*(Amp2(k)-Ampref);       % Controller
            uz(k)          = Cc*xc(1:nc,k) + Dc*(Amp2(k)-Ampref);       % Input to z-piezo [V]
            xcz(1:ncz,k+1) = Acz*xcz(1:ncz,k) + Bcz*(Amp2(k)-Ampref);   % z-piezo response
            z(k)           = Ccz*xcz(1:ncz,k) + Dcz*(Amp2(k)-Ampref);   % Height of scan table [V]
            hz             = h+z(k);                                    % Net height of sample and scan table [V]
        end
        k1 = k1+samplefactor;
        k  = k+1;
    end
    
    %% progress monitoring of simulation
    if i == round(l/100)+prog  
        fprintf(repmat('\b', 1, j)); fprintf('%d%%',round(1/100+prog/l*100));
        prog = i;
        if round(1/100+prog/l*100) == 10
            j = j+1;
        end
    end
end
fprintf(repmat('\b', 1, j)); fprintf('%d%%',100);disp(' ');

%% Sample height estimation for MRM and conventional method
hhat1_miny = height_miny(yc2, uc2, z);                         % Sample height estimation for MRM
hhat1_uz   = height_miny(zeros(size(yc2)), uc2, uz + Ampref);  % Sample height estimation for conventional method

for i = 1:length(uc2)-1                                        % Convert to 2D
hhat_miny(ypos2(i),xpos2(i)) = hhat1_miny(i);
hhat_uz(ypos2(i),xpos2(i))   = hhat1_uz(i);
end
end