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
L = fs / fsx;
x1u = resample(x1, fs, fsx);
x2u = resample(x2, fs, fsx);

t = (0:length(x1u)-1)'/fs;

%% 🎧 MODULACJA DSB-C
y1 = (1 + dA * x1u) .* cos(2*pi*fc1*t);
y2 = (1 + dA * x2u) .* cos(2*pi*fc2*t);
yDSBC = y1 + y2;

%% 🎧 MODULACJA DSB-SC
y1 = dA * x1u .* cos(2*pi*fc1*t);
y2 = dA * x2u .* cos(2*pi*fc2*t);
yDSBSC = y1 + y2;

%% 🎧 MODULACJA SSB-SC
% Filtr Hilberta FIR (np. 101 współczynników)
N = 101;
h = firpm(N-1, [0.05 0.95], [1 1], 'hilbert'); % filtr pasmowo-przepustowy

x1H = filter(h, 1, x1u);
x2H = filter(h, 1, x2u);

% SSB-SC: jedna prawa (USB), jedna lewa (LSB)
y1 = 0.5 * x1u .* cos(2*pi*fc1*t) - 0.5 * x1H .* sin(2*pi*fc1*t); % LSB
y2 = 0.5 * x2u .* cos(2*pi*fc2*t) + 0.5 * x2H .* sin(2*pi*fc2*t); % USB
ySSBSC = y1 + y2;




%% Projektowanie filtru Hilberta (FIR)
N = 201;                   % długość filtru (nieparzysta liczba próbek, typowe dla filtru FIR)
n = -(N-1)/2:(N-1)/2;     % oś próbkowania odpowiedzi impulsowej (indeksy próbek)

% Tworzenie idealnej odpowiedzi impulsowej filtru Hilberta
h_ideal = (1 ./ (pi * n)) .* (1 - cos(pi * n));  % h[n] = 1/(pi*n) * (1 - cos(pi*n))
h_ideal((N+1)/2) = 0;     % ustawienie wartości w n=0 (uniknięcie NaN) zgodnie z definicją

% Zastosowanie okna Hamminga (poprawa charakterystyki filtru)
h = h_ideal .* hamming(N)';  % Okno Hamminga

% Normalizacja filtru Hilberta
h = h / sum(h);  % Normalizacja, aby suma współczynników była równa 1 (zapewnia to stabilność)

%% PRZYGOTOWANIE – DEMODULACJA 
% Projektowanie filtru dolnoprzepustowego (LPF)
lp_cutoff = 4000 / (fs/2); % Przykładowe pasmo: 4 kHz, znormalizowane
lpFilt = fir1(101, lp_cutoff); % Filtr FIR rzędu 101

lowpass_filt = @(sig) filter(lpFilt, 1, sig); % Funkcja filtrująca

% DSB-C
demod_DSB_C1 = resample(lowpass_filt(yDSBC .* cos(2*pi*fc1*t)), fsx, fs);
demod_DSB_C2 = resample(lowpass_filt(yDSBC .* cos(2*pi*fc2*t)), fsx, fs);
demod_DSB_C2 = flipud(demod_DSB_C2);

% DSB-SC
demod_DSB_SC1 = resample(lowpass_filt(yDSBSC .* cos(2*pi*fc1*t)), fsx, fs);
demod_DSB_SC2 = resample(lowpass_filt(yDSBSC .* cos(2*pi*fc2*t)), fsx, fs);
demod_DSB_SC2 = flipud(demod_DSB_SC2);

% SSB-SC
% Tworzymy sygnał analityczny za pomocą filtru Hilberta
hilb = @(x) x + 1i * conv(x, h, 'same');  % sygnał analityczny

% Demodulacja: przesunięcie do pasma podstawowego przez mnożenie z e^{-j2πfct}
sig_anal1 = hilb(ySSBSC) .* exp(-1i * 2 * pi * fc1 * t); % LSB (x1)
sig_anal2 = hilb(ySSBSC) .* exp(-1i * 2 * pi * fc2 * t); % USB (x2)

% Ekstrakcja części rzeczywistej i downsampling
baseband1 = real(sig_anal1);
baseband2 = real(sig_anal2);

demod_SSB_SC1 = resample(baseband1, fsx, fs);
demod_SSB_SC2 = resample(baseband2, fsx, fs);
demod_SSB_SC2 = flipud(demod_SSB_SC2); % dla poprawności kolejności

%% 🎧 ODTWARZANIE
fprintf("\nOdtwarzanie transmisji DSB-C...\n");
disp("Stacja 1 – przed modulacją");
soundsc(x1, fsx); pause();
disp("Stacja 1 – po demodulacji");
soundsc(demod_DSB_C1, fsx); pause();

disp("Stacja 2 – przed modulacją");
soundsc(x2, fsx); pause();
disp("Stacja 2 – po demodulacji");
soundsc(demod_DSB_C2, fsx); pause();

fprintf("\nOdtwarzanie transmisji DSB-SC...\n");
disp("Stacja 1 – przed modulacją");
soundsc(x1, fsx); pause();
disp("Stacja 1 – po demodulacji");
soundsc(demod_DSB_SC1, fsx); pause();

disp("Stacja 2 – przed modulacją");
soundsc(x2, fsx); pause();
disp("Stacja 2 – po demodulacji");
soundsc(demod_DSB_SC2, fsx); pause();

fprintf("\nOdtwarzanie transmisji SSB-SC...\n");
disp("Stacja 1 – przed modulacją");
soundsc(x1, fsx); pause();
disp("Stacja 1 – po demodulacji");
soundsc(demod_SSB_SC1, fsx); pause();

disp("Stacja 2 – przed modulacją");
soundsc(x2, fsx); pause();
disp("Stacja 2 – po demodulacji");
soundsc(demod_SSB_SC2, fsx); pause();
