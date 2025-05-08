clear; close all; clc;

%% Parametry
fs = 400e3; % częstotliwość próbkowania sygnału radiowego
fc1 = 100e3; % nośna 1
fc2 = 110e3; % nośna 2
dA = 0.25; % głębokość modulacji

%% Wczytaj plik audio
[x1, fsx] = audioread('mowa8000.wav');
x2 = flipud(x1); % odwrotnie puszczona mowa

% Normalizacja
x1 = x1 / max(abs(x1));
x2 = x2 / max(abs(x2));

% Nadpróbkowanie
x1u = resample(x1, fs, fsx);
x2u = resample(x2, fs, fsx);
t = (0:length(x1u)-1)'/fs;

%% MODULACJE
% DSB-C
yDSBC = (1 + dA * x1u) .* cos(2*pi*fc1*t) + (1 + dA * x2u) .* cos(2*pi*fc2*t);

% DSB-SC
yDSBSC = dA * x1u .* cos(2*pi*fc1*t) + dA * x2u .* cos(2*pi*fc2*t);

% SSB-SC
% Projektowanie filtru Hilberta
N = 201;
n = -(N-1)/2:(N-1)/2;
h_ideal = (1./(pi*n)) .* (1 - cos(pi*n));
h_ideal((N+1)/2) = 0;
h = h_ideal .* hamming(N)';
hilb = @(x) x + 1i * conv(x, h, 'same');

x1H = filter(h, 1, x1u);
x2H = filter(h, 1, x2u);
ySSBSC = 0.5 * x1u .* cos(2*pi*fc1*t) - 0.5 * x1H .* sin(2*pi*fc1*t) + ...
         0.5 * x2u .* cos(2*pi*fc2*t) + 0.5 * x2H .* sin(2*pi*fc2*t);

%% DEMODULACJA PRZEZ OBWIEDNIĘ (dla każdej modulacji osobno)
% Projektowanie filtrów pasmowoprzepustowych (dla rozdzielenia stacji)
bpFilt1 = designfilt('bandpassfir','FilterOrder',100, ...
    'CutoffFrequency1',95e3,'CutoffFrequency2',105e3,'SampleRate',fs);
bpFilt2 = designfilt('bandpassfir','FilterOrder',100, ...
    'CutoffFrequency1',105e3,'CutoffFrequency2',115e3,'SampleRate',fs);

extract_envelope = @(sig) abs(hilb(sig)); % funkcja obwiedni z Hilberta

% --- DSB-C ---
x_fc1 = filter(bpFilt1, yDSBC);
x_fc2 = filter(bpFilt2, yDSBC);
env1 = extract_envelope(x_fc1);
env2 = extract_envelope(x_fc2);
demod_DSB_C1 = resample(env1, fsx, fs);
demod_DSB_C2 = flipud(resample(env2, fsx, fs));

% --- DSB-SC ---
x_fc1 = filter(bpFilt1, yDSBSC);
x_fc2 = filter(bpFilt2, yDSBSC);
env1 = extract_envelope(x_fc1);
env2 = extract_envelope(x_fc2);
demod_DSB_SC1 = resample(env1, fsx, fs);
demod_DSB_SC2 = flipud(resample(env2, fsx, fs));

% --- SSB-SC ---
x_fc1 = filter(bpFilt1, ySSBSC);
x_fc2 = filter(bpFilt2, ySSBSC);
env1 = extract_envelope(x_fc1);
env2 = extract_envelope(x_fc2);
demod_SSB_SC1 = resample(env1, fsx, fs);
demod_SSB_SC2 = flipud(resample(env2, fsx, fs));

%% ODTWARZANIE
fprintf("\n🔈 Odtwarzanie transmisji DSB-C...\n");
disp("Stacja 1 – przed modulacją"); soundsc(x1, fsx); pause();
disp("Stacja 1 – po demodulacji (obwiednia)"); soundsc(demod_DSB_C1, fsx); pause();

disp("Stacja 2 – przed modulacją"); soundsc(x2, fsx); pause();
disp("Stacja 2 – po demodulacji (obwiednia)"); soundsc(demod_DSB_C2, fsx); pause();

fprintf("\n🔈 Odtwarzanie transmisji DSB-SC...\n");
disp("Stacja 1 – po demodulacji (obwiednia)"); soundsc(demod_DSB_SC1, fsx); pause();
disp("Stacja 2 – po demodulacji (obwiednia)"); soundsc(demod_DSB_SC2, fsx); pause();

fprintf("\n🔈 Odtwarzanie transmisji SSB-SC...\n");
disp("Stacja 1 – po demodulacji (obwiednia)"); soundsc(demod_SSB_SC1, fsx); pause();
disp("Stacja 2 – po demodulacji (obwiednia)"); soundsc(demod_SSB_SC2, fsx); pause();
