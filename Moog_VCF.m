%--------------------------------------------------------------------------
% Ry Currier
% February 24 2022
% 
% Nonlinear Moog VCF 
% Trapezoid FD Scheme with Input Gain, LFO, and Envelope Follower
%--------------------------------------------------------------------------

clc; clear all; close all

%% Input Parameters

SR = 44100;             % Sample Rate (Hz)
f0 = 200;               % Center Cutoff Frequency (Hz)
Tf = 0.2;               % Duration (s)
r = 0.7;                % Feedback Coeff (in [0,1])
gain = 5.0;             % Input Gain (in [1,10])

d = 0.9;                % LFO Depth (proportional to f0)
fLFO = 0.5;             % LFO Frequency (Hz)

de = 4800;              % Envelope Follower Depth (Hz)

maxIter = 3;           % Max Number of Guesses for Newton-Raphson
tol = 1e-10;             % Tolerance for Newton-Raphson

mode = "Audio";         % Toggles Input for 'IR' or 'Audio'

%% Read Input Audio

if mode == "Audio"
    
    [in,SR] = audioread('fingertips.wav');
    u = [0; 0; (1+4*log10(gain))*(sum(in,2)/2)];        % Apply Gain
    Tf = (length(u)-2)/SR;
    env = sqrt(abs(u).^2 + abs(hilbert(u)).^2);         % Get Signal Envelope
    env = env/max(env);                                 % Normalize
    
elseif mode == "IR"
    d = 0;                                              % Make LFO not Happen
    de = 0;                                             % Make Env Follower Not Happen
    env = zeros(Tf*SR+2,1);
end

%% Derived Parameters

k = 1/SR;                               % Time Step (s)
Nf = floor(Tf*SR);                      % Duration (samps)

t = [0:Nf-1]'*k;                        % Time Vector
f = [0:Nf-1]'*SR/Nf;                    % Frequency Vector

fvec = f0*(1 + d*cos(2*pi*fLFO*t));     % Cutoff Frequency Vector (Hz)
om0 = 2*pi*fvec;                        % Angular Cutoff Frequency Vector (rad/s)

%% Vectors and Matrices for DSP Loop

c = [0, 0, 0, 1];                       % Output Reading Vector

x = zeros(4,1);                         % Initialized State Vector
y = zeros(Nf,1);                        % Initialized Output Vector

ga = @(om,r,a1,a2,x1,x2,u0,u1)...
    a1 - x1 + om*k/2*(tanh(a1) + tanh(x1) - tanh(u1 - 4*r*a2) - tanh(u0 - 4*r*x2));
gb = @(om,a1,a2,x1,x2)...
    a2 - x2 + om*k/2*(tanh(a2) + tanh(x2) - tanh(a1) - tanh(x1));

Ja = @(om,a) 1 + om*k/2*(sech(a))^2;
Jb = @(om,a) -om*k/2*(sech(a))^2;
Jc = @(om,r,a,u) 2*r*om*k*(sech(u - 4*r*a))^2;

if mode == "IR"
    u = [0; 1; zeros(Nf-1,1)];          % Input Impulse Response Vector
end

%% DSP Loop

tic;
for n=1:Nf
    
    % Apply Envelope Follower
    
    om0(n) = om0(n) + pi/2*de*(env(n+2) + 2*env(n+1) + env(n));
    
    % Reset Newton-Raphson Parameters
    
    a = x;          % Initial Guess
    iter = 0;       % Iteration Counter
    step = 1;       
    
    while (iter < maxIter) && (norm(step) > tol)
        
        % Evaluate Vector-Valued Function and Jacobian
        
        gmat = [ga(om0(n),r,a(1),a(4),x(1),x(4),u(n),u(n+1));...
            gb(om0(n),a(1),a(2),x(1),x(2));...
            gb(om0(n),a(2),a(3),x(2),x(3));...
            gb(om0(n),a(3),a(4),x(3),x(4))];
        
        Jmat = [Ja(om0(n),a(1)),Jb(om0(n),a(1)),0,0;...
            0,Ja(om0(n),a(2)),Jb(om0(n),a(2)),0;...
            0,0,Ja(om0(n),a(3)),Jb(om0(n),a(3));...
            Jc(om0(n),r,a(4),u(n+1)),0,0,Ja(om0(n),a(4))]';
        
        % Get Step Vector and Update
        
        step = Jmat\gmat;
        a = a - step; 
        iter = iter + 1;
    end

    % Update x
    x = a;

    % Write Sample to Ouput
    y(n) = c*x;
    
end
toc;

if mode == "IR"
    
    % Compute Simulation Transfer Function

    H = abs(fft(y));
    
    % Dunno How to Calculate This Exact Transfer Function

    % Plotting

    loglog(f,H);
    title('Nonlinear Moog VCF Transfer Function')
    subtitle(['r = ' num2str(r) ', f0 = ' num2str(f0)])
    xline(f0,'--')
    xlim([10e0 2*10e3])
    ylim([10e-7 10e1])
    grid on
    xlabel('f (Hz)')
    ylabel('|H| (dB)')
    
elseif mode == "Audio"
    
    % Listen to Ouput
    
    soundsc(y,SR)
    
    % Plotting
    
    subplot(2,1,1)
    %[s1,w1,t1] = 
    spectrogram(u,hamming(128),96,128,SR,'yaxis');
    %imagesc(t1*k^2*Nf,w1*SR^2/Nf,20*log10(abs(s1)),[-30 30]);
    title('Input Sectrogram')
    xlabel('t (s)')
    ylabel('f (Hz)')
    %ax = gca;
    %gca.YScale = 'log';
    colorbar
    caxis([-120 -20])
    
    subplot(2,1,2)
    %[s2,w2,t2] = 
    spectrogram(y,hamming(128),96,128,SR,'yaxis');
    %imagesc(t2*k^2*Nf,w2*SR^2/Nf,20*log10(abs(s2)),[-30 30]);
    title('Output Spectrogram')
    xlabel('t (s)')
    ylabel('f (Hz)')
    %ax = gca;
    %gca.YScale = 'log';
    colorbar
    caxis([-120 -20])
end
