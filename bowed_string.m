%--------------------------------------------------------------------------
% Ry Currier
% 2022-04-12
%
% Bowed Stiff String with Frequency Dependant Loss, Newton-Raphson Solver,
% Interpolated Stereo Output, and Perceptual Controls
%--------------------------------------------------------------------------

clc; clear all; close all

%%%%% flags

plot_on = 0;                % in-loop plotting on (1) or off (0)

%%%%% parameters

% Input Controls

note = 28;                  % MIDI Note Number
stiffness = 0.02;           % Inharmonicity
decay = 5;                  % Fundamental Decay Time (s)
brightness = 3;           % Highest Harmonic Decay Time (s)

N = 100;                    % Total Number of Grid Points

% Bow Parameters

FB = 500;                   % Bow Force / String Mass (m/s^2)
vB = 0.2;                   % Bow Velocity (m/s)
sig = 100;                  % Friction Law Free Parameter (1/m^2)
dur = 1.0;                  % Bowing Duration (s)

% I/O

SR = 44100;                 % sample rate (Hz)
Tf = 5;                     % duration of simulation (s)
tol = 1e-4;                 % Newton-Raphson Tolerance

xi = 0.9;                   % coordinate of excitation (normalised, 0-1)
xoL = 0.2;                  % coordinate of left output (normalized, 0-1)
xoR = 0.25;                 % coordinate of right output (normalised, 0-1)

%% derived parameters

freq = 440*2^((note-69)/12);    % Frequency (Hz)
c = freq*2;                     % wave speed
om = [1 20]*freq;                % T60 frequencies (Hz)
K = sqrt(stiffness)*c/pi;       % Stiffness Constant
T60 = [decay brightness];       % T60 Times

xiom = 1/(2*K^2)*(-c^2 + sqrt(c^4+4*K^2*om.^2));
sig0 = 6*log(10)/(xiom(2)-xiom(1))*(xiom(2)/T60(1) - xiom(1)/T60(2));   % frequency dependent loss parameters
sig1 = 6*log(10)/(xiom(2)-xiom(1))*(1/T60(2) - 1/T60(1));

k = 1/SR;                                                               % time step
Nf = floor(SR*Tf);                                                      % number of time steps

%% Calculate h and N %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hmin = sqrt(0.5*((c^2+4*sig1)*k^2+sqrt((c^4+16*sig1^2+8*sig1*c^2)*k^4+16*K^2*k^2)));    % Minimal Grid Spacing for Stability 
h = hmin+0.0001;                                                                               % Integer Grid Spacing

%% Scheme Matrices and Vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e = ones(N-1,1);                                                        % Ones Vector
Dxx = (1/h^2)*spdiags([e -2*e e], -1:1, N-1,N-1);                       % Sparse 2nd Derivative Matrix

a = 1/(1+sig0*k);                                                       % Coefficient for Neatness
B = a*(2*speye(N-1) + (c^2*k^2+2*sig1*k)*Dxx - K^2*k^2*Dxx*Dxx);        % B Matrix
C = a*((sig0*k-1)*speye(N-1) - sig1*k*Dxx);                             % C Matrix

xi_int = floor(xi * N);                                                 % Integer Input Grid Index
xoL_int = floor(xoL * N);                                               % Integer Left Output Grid Index
xoL_frac = xoL*N - xoL_int;                                             % Fractional Left Output Grid Index
xoR_int = floor(xoR * N);                                               % Integer Left Output Grid Index
xoR_frac = xoR*N - xoR_int;                                             % Fractional Left Output Grid Index

J = k^2/h*[zeros(xi_int-1,1); 1; zeros(N-xi_int-1,1)];                  % J Vector
J = sparse(J);                                                          % Sparse J Vector
cvL = [zeros(1,xoL_int-1), 1-xoL_frac, xoL_frac, zeros(1,N-xoL_int-2)]; % Left Channel c Vector
cvR = [zeros(1,xoR_int-1), 1-xoR_frac, xoR_frac, zeros(1,N-xoR_int-2)]; % Left Channel c Vector
cv = sparse([cvL; cvR]);                                                % Sparse c Vector

dur_int = floor(dur*SR);                                                % Bowing Duration (samples)
FBv = FB*[ones(dur_int,1); zeros(Nf-dur_int,1)];                        % Bow Force Vector

%% initialise scheme variables

u2 = zeros(N-1,1);              % state
u1 = u2;                        % state
u = u2;                         % state
ulast = u2;                     % NR Iterative state

y = zeros(Nf,2);                % output
xax = (1:N-1)';               % x-axis for plotting

%% main loop

tic
for n=1:Nf
    
    b = B*u1 + C*u2;            % Undriven Scheme Update
    eps = 1;                    % Initial Error for NR
    
    % Newton-Raphson
    
    while eps>tol
        
        % Coefficients for Neatness
        
        eta = 1/(2*k)*(ulast(xi_int) - u2(xi_int)) - vB;
        phi = J*FBv(n)*exp(0.5)*sqrt(2*sig)*eta*exp(-sig*eta^2);
        
        % NR Update
        
        u = ulast - (ulast - b - phi)./(1 - phi + sig*eta^2*phi);
        
        % Compute Error and Shift State
        
        eps = abs(u-ulast);
        ulast = u;
        
    end
    
    % read output
    
    y(n,:) = cv*u;
    
    % plot
    
    if (plot_on==1)
        % draw current state
        figure(1)
        plot(xax, u, 'k');
        xlabel('x (m)')
        ylabel('u (m)')
        axis([0 N -0.005 0.005])
        drawnow
    end
    
    % shift state
    uvec(n) = u2(xi_int); 
    u2 = u1;
    u1 = u;
    
end
toc

%% play sound

soundsc(y,SR);

%% plot spectrum

figure(2)

yfft = 10*log10(abs(fft(y)));
plot([0:Nf-1]*Tf/Nf, sum(y,2), 'k')
xlabel('Time (s)')
ylabel('Magnitude')
title('Output')
