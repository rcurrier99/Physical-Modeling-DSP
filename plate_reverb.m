%--------------------------------------------------------------------------
% Ry Currier
% March 22 2022
%
% Plate Reverb 
% Pruned Modal Model with Material Selection, Frequency-Dependent Loss,
% Size, Tension, Thickness, Input/Output Location, and Wet/Dry Controls
%--------------------------------------------------------------------------

clc; clear all; close all

%% Read Audio

[in, SR] = audioread("fingertips.wav");
in = sum(in,2);

%% Parameters

T = 700;                                    % Tension (N)
H = 5e-4;                                   % Thickness (m)
E = [9.6e9, 6.8e10, 1.1e11, 2.02e11];       % Young's modulus (Pa) (wood, aluminum, copper, steel)
rho = [700, 2600, 8960, 7860];              % Density (kg/m^3) (wood, aluminum, copper, steel)
v = [0.4, 0.32, 0.33, 0.29];                % Poisson's Ratio  (wood, aluminum, copper, steel)
Lx = 2;                                     % X Length (m)
Ly = 1;                                     % Y Length (m)
T60 = [5 0.5];                              % T60 times (s)

mat = 0.5;                                  % Material Selector (0-1)
wetdry = 0.5;                               % Wet/Dry Control (0-1)

%% Material Interpolation

mat = 3*mat+1;                              % Material Selector in Range 1-4
P = zeros(4,1);                             % Store Interpolation Coefficients

% Compute Interpolation Coefficients
for j=1:4
    l = [1:4];
    l(j) = [];
    P(j) = prod(mat-l)/prod(j-l);
end

% Interpolate Physical Quantities
E = E * P;
rho = rho * P;
v = v * P;

%% Error Checking

if T < 0 | H < 0 | E < 0 | rho < 0 | v < 0 | T60(1) < 0 | T60(2) < 0 | Lx < 0 | Ly < 0
    error("Parameters Cannot Be Negative")
end

%% Derived Parameters

k = 1/SR;                               % Time Step
K = sqrt(E*H^2/(12*rho*(1-v^2)));       % Stiffness Parameter
c = sqrt(T/(rho*H));                    % Wave Speed

z = T60(1)*SR;                          % Fundamental T60 Time (samps)
in = [in; zeros(z,1)];                  % Backpad Input by T60

Nf = length(in);                        % Simulation Length (samps)
Tf = Nf/SR;                             % Simulation Time (s)

%% Maximum Modal Frequency

w_max = 2/k;                                                % Maximum Frequency
beta_max = sqrt((-c^2+sqrt(c^4+4*K^2*w_max^2))/(2*K^2));    % Maximum Wavenumber

Mx = floor(Lx/pi * sqrt(beta_max^2-(pi/Ly)^2));             % Maximum x Index
My = floor(Ly/pi * sqrt(beta_max^2-(pi/Lx)^2));             % Maximum y Index
M = Mx*My;                                                  % Maximum Modal Pair Index

%% Modal Index Pairs

[mx,my] = meshgrid([1:Mx],[1:My]);                          % Create Pair of Index Matricies
mx = reshape(mx,[Mx*My,1]); my = reshape(my,[Mx*My,1]);     % Reshape to Pair of Column Vectors

beta = sqrt((mx*pi./Lx).^2 + (my*pi./Ly).^2);               % Wavenumber Vector
omega = sqrt(c^2*beta.^2 + K^2*beta.^4);                    % Undamped Modal Frequencies

%% Pruning Round 1
 
% Round Frequencies to 3 Decimal Places, Keep Unique Elements
[omega, iomega1, iomega2] = unique(omega-mod(omega,0.001),'stable');

% Keep Corresponding Entries of Beta, mx, my
beta = beta(iomega1);
mx = mx(iomega1);
my = my(iomega1);

stab = beta < beta_max;                                     % 1 if stable, 0 else

%% Loss

sig_coeff = 6*log(10)/(beta(end)^2-beta(1)^2);              % Coefficient for Brevity
sig0 = sig_coeff*(beta(end)^2/T60(1)-beta(1)^2/T60(2));     % Sigma0
sig1 = sig_coeff*(1/T60(2)-1/T60(1));                       % Sigma1
sigma = sig0 + sig1*beta.^2;                                % Loss Vector

if size(nonzeros(sigma > 0),1) ~= size(sigma,1)
    error("Loss Vector Contains Non-positive Elements!")
end

%% I/O

xi = [0.1,0.1];                                             % Input Location (Normalized 0-1)
xo1 = [0.8,0.9];                                            % Output Location 1 (Normalized 0-1)
xo2 = [0.9,0.8];                                            % Output Location 2 (Normalized 0-1)

IOc = 2/sqrt(Lx*Ly);                                        % Coefficient for Brevity
phii = IOc*sin(mx*pi*xi(1)/Lx).*sin(my*pi*xi(2)/Ly);        % Mode Shapes at Input
phio1 = IOc*sin(mx*pi*xo1(1)/Lx).*sin(my*pi*xo1(2)/Ly);     % Mode Shapes at Output 1
phio2 = IOc*sin(mx*pi*xo2(1)/Lx).*sin(my*pi*xo2(2)/Ly);     % Mode Shapes at Output 2

izero = phii > 1e-13;                                       % Modes Excited at Input
ozero1 = phio1 > 1e-13;                                     % Modes Contributing to Output 1
ozero2 = phio2 > 1e-13;                                     % Modes Contributing to Output 2
zero = izero.*or(ozero1,ozero2);                            % Modes Excited and Contributing to Either Output

phii = nonzeros(zero.*stab.*phii);                          % Stable Nonzero Mode Shapes at Input
phio1 = nonzeros(zero.*stab.*phio1);                        % Stable Nonzero Mode Shapes at Output 1
phio2 = nonzeros(zero.*stab.*phio2);                        % Stable Nonzero Mode Shapes at Output 1

out = zeros(Nf,2);                                          % Initial Output

%% Pruning Round 2

% Wavenumbers Frequencies and Losses Corresponding to Stable Nonzero Modes
beta = nonzeros(zero.*stab.*beta);
omega = nonzeros(zero.*stab.*omega);
sigma = nonzeros(zero.*stab.*sigma);

%% Scheme Variables

m = length(beta);                                           % Number of Modes
p2 = zeros(m,1);                                            % State
p1 = p2;
    
B = (2-k^2*omega.^2)./(1+k*sigma);                          % Update Coeffs
C = (k*sigma-1)./(1+k*sigma);

%% Main Loop

tic
for n=1:Nf
    
    % Update State
    p = B.*p1 + C.*p2 + phii*in(n);                            
    
    % Shift State
    p2 = p1;    
    p1 = p;
    
    % Add to Output
    out(n,1) = phio1'*p;
    out(n,2) = phio2'*p;
    
end
toc

%% Listen to Output

wet = diff(out./max(out));
wet = [[0; wet(:,1)], [0; wet(:,2)]];
wet = wet./max(wet);

dry = [in,in];
dry = dry./max(dry);

audio = ((1-wetdry)*wet + wetdry*dry)/2;
soundsc(audio,SR);

%% Plotting

fax = [0:Nf-1]*SR/Nf;
xax = fax*Tf/SR;
audiofft = 20*log10(abs(fft(sum(audio,2))));

subplot(2,1,1)
plot(xax,sum(audio,2))
xlim([0 Tf])
xlabel('Time (s)')
ylabel('Magnitude')
title('Output')

subplot(2,1,2)
semilogx(fax,audiofft);
xlim([1 SR/2]);
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Transform of Output')
