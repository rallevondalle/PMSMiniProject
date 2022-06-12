%2D Wave Guide
%clear;
%% Initialise variables
fs = 44100;             % sample rate [Hz]
k = 1 / fs;             % time step [s]
M = 1;                  % mass
K = 1;
kappa = 0.01;           % set to smaller for spring, like kappa = 0,006
c = 300;                % speed of sound in air, 342 m/s
Lx = 2;                 % length of x axis
Ly = Lx;                % length of y axis
sigma0 = 0.8;           %
sigma1 = 0.005;         %
h = sqrt(2*(c.^2*k.^2 + 4*sigma1*k)); %Initialize grid space, calculated from the Courant-number (Stability condition)
Nx = floor(Lx/h);       %
Ny = floor(Ly/h);       %
h = Lx / Nx;            % assuming grid space is equal for both dimensions

lengthSound = 5*fs;    % length of the simulation in seconds

f0x = c / (2*Lx);       % fundamental frequency [Hz]
lambda = (c*k)/h;       % CFL condition, Courant-Friedrichs-Lewy
muSq = ((kappa*k)/h^2)^2; %stability constant, 0:1
pickupPosition = floor(0.2*(fs/Nx)); %Samplerate dependet pickup position, close to left (N = 0) position
pickSize = 1;           %Size of pick from 1-10 (smaller number yields higher frequencies and vice versa)

%% Stability check
stability = 1;
% if h >= c*k
%     stability = 1;
% else
%     stability = 0;
% end

%omega0 = 2 * pi * f0; % angular (fundamental) frequency [Hz]
%M = 100; % mass [kg]
%K = omega0^2 * M; % spring constant [N/m]

%Set boundary conditions (clamped, simply supported, free [wind instruments])
%Clamped l = {2,...,N-2}
%Unexpanded FD scheme is deltatt*uln = -k^2deltaxxxx*unl
%Expanded dxxxxunl= 1/h^4(ul+2 - 4ul+1 + 6ul - 4ul-1 + ul-2)
%Boundary condition uln = dxxuln = 0 at l=0,N --> un-1 = -un1 and uN+1 =
%-uN-1


% initialise matrices, not vectors in this case since we are in two
% dimension
uNext   = zeros(Nx+1,Ny+1);
u       = zeros(Nx+1,Ny+1);
uPrev   = zeros(Nx+1,Ny+1);
range   = zeros(Nx+1,Ny+1); %total grid space with boundary conditions
out     = zeros(lengthSound,1);


% plucking string, initial conditions (u0 = 1, d/dt u0 = 0)
% pickSizeScaled = floor(pickSize*0.3*(fs/Nx)) %Pick size in samples dependent on fs and N
% u(pickSizeScaled:2*pickSizeScaled) = hann(pickSizeScaled+1); %Raised cosine with hann function
% uPrev = u;

%Initial conditions, vector lenght of x, add raised cosine to middle,
%transpose and vector multiply to get a sqaure matrix size of your plate
pluckVector = zeros(Nx+1,1);
pluckPosition = floor((Nx/2)+5);
pluckSize = 5;
pluckVector(pluckPosition-pluckSize:pluckPosition+pluckSize) = rkHannWin(2*pluckSize);
pluckMatrix = pluckVector * transpose(pluckVector);
u = pluckMatrix;
uPrev = u;
%figure;
%mesh(u);

%[X,Y] = meshgrid(-1:0.02:1);
X = linspace(-1,1,Nx+1)'; Y = X;
%Z = uNext;
%mesh(X,Y,Z)

visualizeOn = 0;        %visualize plate behaviour

% Simulation loop
if stability == 1
    %disp("System stable");
    %disp(pickupPosition);
    for n = 3:lengthSound %TIME ITERATION, LAPLACIAN DISCRETISATION

        uNext(2:Nx,2:Ny) = -uPrev(2:Nx,2:Ny) + 2*u(2:Nx,2:Ny) + lambda*lambda*(u(3:Nx+1,2:Ny) + ...
                           u(2:Nx,3:Ny+1) - 4*u(2:Nx,2:Ny) + u(1:Nx-1,2:Ny) + u(2:Nx,1:Ny-1)) + ...
                           ((2*sigma1*k)/(h.^2)) * (u(3:Nx+1,2:Ny) + u(1:Nx-1,2:Ny) + u(2:Nx,3:Ny+1) + ...
                           u(2:Nx,1:Ny-1) - 4*u(2:Nx,2:Ny) - uPrev(3:Nx+1,2:Ny) - uPrev(1:Nx-1,2:Ny) - ...
                           uPrev(2:Nx,3:Ny+1) - uPrev(2:Nx,1:Ny-1) + 4*uPrev(2:Nx,2:Ny)) ; %with frequency dependent damping
    
        if visualizeOn == 1
            Z = uNext;
            mesh(X,Y,Z);
            zlim([-0.5 0.5]);
            drawnow;
        end
        
        % Update output at specific position in space, simulating a pick-up (single point) 
        out(n,:) = u(pluckPosition-10,pluckPosition-10);

        % Kinetic energy
        %kinEnergy(n) = M / 2 * (1/k * (u - uPrev))^2;
        
        % Potential energy
        %potEnergy(n) = K / 2 * u * uPrev;
        
        % Total energy (Hamiltonian)
        %totEnergy(n) = kinEnergy(n) + potEnergy(n);
            
        % Update system states    
        uPrev = u;
        u = uNext;    
    end
else
    disp("System unstable, check variables");
end

%% OUTPUT TO WAVEFILE
mkdir(sprintf('audioExports'));

% TIME STAMP FOR FILE OUTPUT
filename = sprintf('audioExports/plate_%ihz_%s.wav',fs,datestr(now,'yyyymmdd-HHMMSS'));
disp(['Exported to ./' filename ' at ' num2str(fs) 'Hz']);

%%audiowrite(filename,out,fs);

%% PLOT

figure("Name","Time domain plot");
plot(out);
soundsc(out,fs); 

%% PLOT SPECTROGRAM

wlen  = 1024;                   % window length (recomended to be power of 2)
hop   = wlen/4;                 % hop size (recomended to be power of 2)
nfft  = 4096;                   % number of fft points (recomended to be power of 2)
N     = length(out);            % N data points
tWave = (0:N-1)/(fs);           % time vector for waveform plot

% perform STFT
win = blackman(wlen);
[S, f, t] = stft(out, win, hop, nfft, fs);

% calculate the coherent amplification of the window
C = sum(win)/wlen;

% take the amplitude of fft(x) and scale it, so not to be a
% function of the length of the window and its coherent amplification
S = abs(S)/wlen/C;

% correction of the DC & Nyquist component
if rem(nfft, 2)                     % odd nfft excludes Nyquist point
    S(2:end, :) = S(2:end, :).*2;
else                                % even nfft includes Nyquist point
    S(2:end-1, :) = S(2:end-1, :).*2;
end

% convert amplitude spectrum to dB (min = -120 dB)
S = 20*log10(S + 1e-6);

% plot the spectrogram
figure(3)
subplot(2,1,1);
surf(t, f, S)
shading interp
axis tight
view(0, 90)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Time, s')
ylabel('Frequency, Hz')
title('Spectrogram')

hcol = colorbar;
set(hcol, 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel(hcol, 'Magnitude, dB')

subplot(2,1,2);
plot(tWave,out,'black');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Time, s')
ylabel('Amplitude')
title('Amplitude Waveform')
ylim([-0.2 0.2]);

%spectrogram(out,512,64,512, fs, 'yaxis');
%plot(out_fft);


% Window functions
function [window] = rkHammingWin(size);
if size > 0
    M = size;
    w = .54 - .46*cos(2*pi*(0:M)'/(M-1));       %for 1 based systems 0:M, for 0 based systems 0:M-1
else
end
end

function [window] = rkHannWin(size)
if size > 0    
    M = size;    
    window = 0.5 * (1 - cos(2*pi*(0:M)'/M));    %for 1 based systems 0:M, for 0 based systems 0:M-1
else
end
end


% STFT FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Short-Time Fourier Transform            %
%               with MATLAB Implementation             %
%                                                      %
% Author: Ph.D. Eng. Hristo Zhivomirov        12/21/13 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [STFT, f, t] = stft(x, win, hop, nfft, fs)

% function: [STFT, f, t] = stft(x, win, hop, nfft, fs)
%
% Input:
% x         - signal in the time domain
% win       - analysis window function
% hop       - hop size
% nfft      - number of FFT points
% fs        - sampling frequency, Hz
%
% Output:
% STFT      - STFT-matrix (only unique points, time 
%             across columns, frequency across rows)
% f         - frequency vector, Hz
% t         - time vector, s

% representation of the signal as column-vector
x = x(:);

% determination of the signal length 
xlen = length(x);

% determination of the window length
wlen = length(win);

% stft matrix size estimation and preallocation
NUP = ceil((1+nfft)/2);     % calculate the number of unique fft points
L = 1+fix((xlen-wlen)/hop); % calculate the number of signal frames
STFT = zeros(NUP, L);       % preallocate the stft matrix

% STFT (via time-localized FFT)
for l = 0:L-1
    % windowing
    xw = x(1+l*hop : wlen+l*hop).*win;
    
    % FFT
    X = fft(xw, nfft);
    
    % update of the stft matrix
    STFT(:, 1+l) = X(1:NUP);
end

% calculation of the time and frequency vectors
t = (wlen/2:hop:wlen/2+(L-1)*hop)/fs;
f = (0:NUP-1)*fs/nfft;

end
