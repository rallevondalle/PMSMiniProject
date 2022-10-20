%2D Wave Guide
%clear;
%% Initialise variables
fs = 44100;             % sample rate [Hz]
k = 1 / fs;             % time step [s]
M = 1;                  % mass
K = 1;                  % spring constant
T = 1000000;            % tenstion per unit length T (in N / m), 'spring constant'
rho = 7850;             % density in kg/m^3
H = 0.0005;             % thickness in meters
c = sqrt(T / (rho*H));  % speed of wave in medium
Lx = 2;                 % length of x axis in meters
Ly = Lx;                % length of y axis in meters

kappa = 0.006;          % set to smaller for spring, like kappa = 0,006 [curvature]
sigma0 = 1.0;           % frequency independent damping
sigma1 = 0.005;         % frequency dependent damping (set very low for spring)

h = sqrt(2*(c.^2*k.^2 + 4*sigma1*k)); %Initialize grid space, calculated from the Courant-number (Stability condition)
Nx = floor(Lx/h);       % number of grid points x-axis (l)
Ny = floor(Ly/h);       % number of grid points y-axis (m)
h = min(Lx / Nx);       % assuming grid space is equal for both dimensions

f0x = c / (2*Lx);       % fundamental frequency [Hz] ONLY FOR 1D STRING
lambda = (c*k)/h;       % CFL condition, Courant-Friedrichs-Lewy, lambda lessequals 1/sqrt(2)
muSq = ((kappa*k)/h^2)^2;% stability constant, 0:1
pickupPosition = floor(0.2*(fs/Nx)); %Samplerate dependet pickup position, close to left (N = 0) position
pluckPositionSkew = 25; % move pluck position away from center of plate [5 45]


objectName = '2Dplate'; % used for audio and figure outputs
analysisLength = 10;     % number of seconds to compute
lengthSound = analysisLength*fs;% length of the simulation in seconds
pickSize = 15;          % size of pick from 1-50 (smaller number yields higher frequencies and vice versa)

% Stability check
if h >= sqrt(2)*c*k
    stability = 1;
    fprintf('System stable: h = %f >= c*k = %f\n',h,c*k);
    fprintf('Proceeding with calculations\n')
else
    stability = 0;
    msg = sprintf('System unstable! h = %f >= c*k = %f\n',h,c*k');
    error(msg);
end

%omega0 = 2 * pi * f0; % angular (fundamental) frequency [Hz]
%M = 100; % mass [kg]
%K = omega0^2 * M; % spring constant [N/m]

% Set boundary conditions (clamped, simply supported, free [wind instruments])
%Clamped l = {2,...,N-2}
%Unexpanded FD scheme is deltatt*uln = -k^2deltaxxxx*unl
%Expanded dxxxxunl= 1/h^4(ul+2 - 4ul+1 + 6ul - 4ul-1 + ul-2)
%Boundary condition uln = dxxuln = 0 at l=0,N --> un-1 = -un1 and uN+1 =
%-uN-1


% Initialise 2D matrices
uNext   = zeros(Nx+1,Ny+1);
u       = zeros(Nx+1,Ny+1);
uPrev   = zeros(Nx+1,Ny+1);
range   = zeros(Nx+1,Ny+1); %total grid space with boundary conditions
out     = zeros(lengthSound,1);


% Excite system
% Initial conditions, vector lenght of x, add raised cosine to middle,
% transpose and vector multiply to get a sqaure matrix size of your plate
pluckVector = zeros(Nx+1,1);
pluckPosition = floor((Nx/2)+pluckPositionSkew); % pluckPositionSkew moves the excitation away from center
pluckSize = pickSize;
pluckVector(pluckPosition-pluckSize:pluckPosition+pluckSize) = rkHannWin(2*pluckSize);
pluckMatrix = pluckVector * transpose(pluckVector);
u = pluckMatrix;
uPrev = u;

% Plot plucking matrix
figure(1);
mesh(u);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Nx')
ylabel('Ny')
zlabel('Amplitude')
title(sprintf('Pluck matrix, pick size = %i',pickSize))

%%

visualizeOn = 1;        % visualize plate behaviour

% Simulation loop
X = linspace(-1,1,Nx+1)'; Y = X;

if stability == 1
    %fprintf('Pickup position = %i\n',pickupPosition);
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
    disp("System unstable, check variables!");
end

soundsc(out,fs);

% PLOT WAVEFORM
figure(2);
plot(out);
title('Time domain plot')
xlabel('Time (s)')
ylabel('Amplitude (dB)')
%% OUTPUT TO WAVEFILE
myFolder = 'audioExports';
if not(isfolder(myFolder))
    sprintf('Created directory %s',myFolder)
    mkdir(myFolder)
else
    sprintf('Directory "%s" already exists',myFolder)
end

% TIME STAMP FOR FILE OUTPUT
filename = sprintf('%s/%s_%ihz_%s.wav',myFolder,objectName,fs,datestr(now,'yyyymmdd-HHMMSS'));
disp(['Exported to ./' filename ' at ' num2str(fs) 'Hz']);

audiowrite(filename,out,fs);   % uncomment to print output matrix to wavefile

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
ylim([-1 1]);


%Print spectrogram to file
myFolder = 'figures';
if not(isfolder(myFolder))
    sprintf('Created directory %s',myFolder)
    mkdir(myFolder)
else
    sprintf('Directory "%s" already exists',myFolder)
end

filename = sprintf('%s/%s_spectrogram_%ihz_%isec',myFolder,objectName,fs,analysisLength);
print(filename,'-dpng');
disp(['Exported to ./' filename ' at ' num2str(fs) ' Hz']);


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
