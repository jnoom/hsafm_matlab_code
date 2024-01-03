clearvars
rng(2022)
addpath('functions')
addpath('FISTA-master/FISTA-master/');
addpath('FISTA-master/FISTA-master/utils/');

%% Construct sample height h_sample and the height considering the tip radius h_radius
lines     = 20;       % lines per frame
scanspeed = 200*10^3; % scan speed [nm/s]
scansize  = 24;       % scan size [nm]
f_drive   = 4*10^5;   % drive frequency of the cantilever [Hz]
steps     = 8;        % max. number of steps in sample (Note: make sure that "f_drive/scanspeed*scansize/steps/2" is an integer)
c         = 40;       % hardness of sample (large value can induce instability)
height    = 1;        % step height [nm]
base      = -2;       % base height [nm]
[h_sample,h_radius] = construct_sample(base,height,f_drive,scanspeed,scansize,steps,lines);  % sample height [nm]

%% Controller
Kp     = 0.1;       % proportional gain
Ki     = 7.5*10^4;  % integral gain
Kd     = 0;         % derivative gain
alphak = 5*10^5;    % derivative cutoff
ts     = 10^-7;     % sampling time [s]
controller = ss(c2d(tf([Kp Ki],[1 0])+tf([Kd 0],[1 alphak]),ts));  % Controller

%% z-piezo
f_zpiezo    = 4*10^4;  % resonance frequency of z-piezo [Hz]
beta_zpiezo = 0.5;     % damping factor of z-piezo
zpiezo      = ss(c2d(tf((f_zpiezo*2*pi)^2,[1 beta_zpiezo*(f_zpiezo*2*pi) (f_zpiezo*2*pi)^2]),ts));  % z-piezo dynamics

%% Simulate cantilever response
nm2V      = 0.08/2;      % nanometer to Volts conversion
amp_drive = 0.48;        % drive amplitude of the cantilever [V]
Ampref    = 1.6*nm2V;    % reference amplitude for feedback [V]
phase     = -pi/3;       % phase delay for input [rad]
fcut      = f_drive*0.4; % cutoff frequency of low-pass filter in LIA

load('data/identification','sys','Rnu','Reta','S','tau_opt');  % load identified cantilever dynamics
sys = c2d(d2c(sys),ts);

uc = amp_drive*sin(f_drive*2*pi*[0:ts:scansize/scanspeed*lines*2] - phase)';  % input to cantilever
x0 = zeros(length(sys.B),1);
[yc,xpos,ypos,z,hhat_miny,hhat_uz] = ...
    cantilever_response(sys,x0,uc,Reta,Rnu,h_radius*nm2V,c,scansize,scanspeed,lines,f_drive,fcut,controller,zpiezo,Ampref);

%% STATESON
Q            = 0.10^2;  % estimated variance of pulses
T            = 300;     % nreak the problem down into smaller parts
lambdafactor = 1;       % factor for tuning lambda
epsfactor    = 1;       % factor for epsilon
err          = 0.01;    % minimum magnitude of pulses

tic;
[Fhhat,yhat] = stateson(uc,yc,sys,length(sys.B),Rnu,Reta,Q,T,lambdafactor,epsfactor,err);  % estimate t/s-interaction and cantilever deflection
Tcomp        = toc;

hhat1_STAT = height_miny(yhat, uc, z);  % determine height from estimated cantilever deflection
for i = 1:length(uc)-1
    hhat_STAT(ypos(i),xpos(i)) = hhat1_STAT(i); % convert to 2D
end

%% Benchmark (Kalman Filter Formulation)
Td         = 6*10^-6;  % decay of impulse response in seconds
Rdelta     = 30;       % variance of delta
Rgamma     = Rnu*50;   % variance of gamma
n          = 5;        % order of system with state theta; at least 2
pulse_time = 8;        % time delay for pulse occurrences

tic
Fhhat_KFF = KFF(yc, uc, sys, Td, Rdelta, Rgamma, n, pulse_time);  % estimate t/s-interaction
Tcomp_KFF = toc;

yhat_KFF = lsim(sys,uc+Fhhat_KFF,0:ts:(length(uc)-1)*ts);         % estimate cantilever deflection

hhat1_KFF = height_miny(yhat_KFF, uc, z);                         % determine height from estimated cantilever deflection
for i = 1:length(uc)-1
    hhat_KFF(ypos(i),xpos(i)) = hhat1_KFF(i);                     % convert to 2D
end

%% Plots
figure;subplot(2,3,1);imagesc([0 scansize],[0 scansize],h_sample);
xlabel({'x [nm]';'Ground truth'});ylabel('y [nm]');
subplot(2,3,2);imagesc([0 scansize],[0 scansize],hhat_uz)
xlabel({'x [nm]';'Conventional'});ylabel('y [nm]')
subplot(2,3,4);imagesc([0 scansize],[0 scansize],hhat_KFF)
xlabel({'x [nm]';'KFF'});ylabel('y [nm]');
subplot(2,3,5);imagesc([0 scansize],[0 scansize],hhat_STAT);
xlabel({'x [nm]';'SBR'});ylabel('y [nm]')
subplot(2,3,6);imagesc([0 scansize],[0 scansize],hhat_miny)
xlabel({'x [nm]';'MRM'});ylabel('y [nm]');