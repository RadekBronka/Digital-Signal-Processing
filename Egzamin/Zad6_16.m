% Wczytanie sygnału
load('ECG100.mat');  % lub 'ECG100 (1).mat'
x = double(val(1, :));  % pierwszy kanał jako double

fs = 360;               % częstotliwość próbkowania w Hz
t = (0:length(x)-1)/fs;

% Centrowanie sygnału
x = x - mean(x);

% FFT widmo
N = length(x);
X = abs(fft(x)); %liczymy transforamtę fouriera
f = (0:N-1)*(fs/N);

figure;
subplot(2,1,1);
plot(f(), X());
xlim([0 10]);%zakres
title('Widmo FFT sygnału EKG');
xlabel('Częstotliwość [Hz]');
ylabel('|X(f)|');
grid on;


% --- Autokorelacja z użyciem FFT ---
R = ifft(abs(fft(x)).^2, 'symmetric'); %obliczamy widmo mocy i robimy ifft
R = R / max(R);  % normalizacja
lags = (0:N-1)/fs; %wektor opóźnień

%szukamy
min_hr = 40; max_hr = 180; %ograniczenia fizjologiczne
min_lag = round(fs * 60 / max_hr);  % ~0.33 s
max_lag = round(fs * 60 / min_hr);  % ~1.5 s

[~, peak_idx] = max(R(min_lag:max_lag));
true_idx = peak_idx + min_lag - 1;
period = lags(true_idx);
bpm = 60 / period;

subplot(2,1,2);
plot(lags(1:round(3*fs)), R(1:round(3*fs)));
xline(period, '--r', sprintf('T = %.2f s', period));
title('Autokorelacja (FFT)');
xlabel('Opóźnienie [s]');
ylabel('Znormalizowana wartość');
grid on;

fprintf('Szacowana liczba uderzeń na minutę: %.2f BPM\n', bpm);
