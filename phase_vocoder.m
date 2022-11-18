%--------------------------------------------------------------------------
% Ry Currier
%
% November 8 2021
%
% Phase Vocoder Pitch-Shifter/Time-Stretcher
%--------------------------------------------------------------------------

clc; clear all; close all

%% Input Controls

PL = 1;                                 % Phase-Locking on/off
PS = 100;                               % Pitch-Shift Amount (Hz) 
Q = 5;                                  % Time-Stretch Factor

T = 25;                                 % Window Length (ms)
O = 0.95;                               % Overlap Factor

file = 'fingertips.wav';                % Audio File to Process

%% Read Audio File

[x,Fs] = audioread(file);
x = 0.5*sum(x,2);

%% Derived Constants

N = floor(Fs/T);                        % Frame Length (samples)
NFFT = 2^ceil(log2(N));                 % FFT Length (samples)
HA = floor(N*(1-O));                    % Analysis Hop Size (samples)
HS = floor(Q*HA);                       % Synthesis Hop Size (samples)

%% Zero Pad Signal

x = [zeros(1,N) x'];                    % Front Pad 1 Frame of Zeros 
L = length(x);                          % Signal Length (samples)
NF = ceil(L/HA);                        % Number of Frames
z = NF*HA-L; 
x = [x zeros(1,z)];                     % 

%% Generate Window

t = [0:N-1];                            % Time Vector
win = (1-cos(2*pi*t/N));                % Hann Window

%% Synthesis

y = zeros(1,NF*HS);                     % Output Vector
phim = zeros(1,NFFT/2+1);               % Phase Vectors
phim1 = zeros(1,NFFT/2+1);
Ym = zeros(1,NFFT/2+1);                 % Transform Vectors
Ym1 = zeros(1,NFFT/2+1);
thetam = zeros(1,NFFT/2+1);
thetam1 = zeros(1,NFFT/2+1);
omegak = [0:NFFT/2]/(NFFT/2+1)*(2*pi);  % Frequency Vector
q = max([N/HA N/HS]);                   % Max Loop Length

%% Loop Over Frames

for m=1:NF-q
    
    X = fft(0.5*exp(i*PS*t/N).*x(m*HA+1:m*HA+N).*win,NFFT);       % Windowed FFT of Current Frame 
    
    Xmag = abs(X);      Xmag = Xmag(1:NFFT/2+1);
    Xang = angle(X);    phim1 = Xang(1:NFFT/2+1);
    
    Omegam1 = omegak + ppa(phim1-phim-omegak*HA)/HA;            % Frequency Estimation
    thetam1 = (1-PL)*thetam + HS*Omegam1;                       % Phase Estimation
    
    Zm1 = PL*angle(Ym(1:NFFT/2+1) - [Ym(2:NFFT/2+1) 0] - [0 Ym(1:NFFT/2)]);
    Ym1 = Xmag.*exp(i*(thetam1 + Zm1)); 
    Ym1 = [Ym1(1:NFFT/2+1) conj(flip(Ym1(2:NFFT/2)))];          % Hermitian Symmetry 
    ym1 = ifft(Ym1);
    
    y(m*HS+1:m*HS+N) = y(m*HS+1:m*HS+N) + ym1(1:N).*win;        % Add Frame to Output
    
    phim = phim1;                                               % Recurse
    thetam = thetam1; 
    Ym = Ym1;

end

soundsc(real(y),Fs)                                             % Hear Output

%% Spectrogram

subplot(2,1,1)
spectrogram(x,win,floor(O*N),N,Fs,'yaxis')
subtitle('Original Signal')
xlabel('Time (s)')
ylabel('Frequency (kHz)')
ax = gca;
ax.YLim=[1 20];
%ax.YScale='log';

subplot(2,1,2)
spectrogram(real(y),win,floor(O*N),N,Fs,'yaxis')
subtitle('Time-Stretched Signal')
xlabel('Time (s)')
ylabel('Frequency (kHz)')
ylim([1 20])
ax = gca;
ax.YLim=[1 20];
%ax.YScale='log';

%% Principal Phase Angle Function
% Eats Phase Angle, Returns Principal Branch Phase Angle
function phi = ppa(x)
    phi = mod(pi+x,2*pi) - pi;
end