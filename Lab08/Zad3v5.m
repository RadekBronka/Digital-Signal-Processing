clear all; close all; clc;

%% 1. Parametry sygnałów
fs_target = 48000; % Docelowa częstotliwość próbkowania
t = 1; % Czas trwania sygnału

% Częstotliwości i próbkowania sygnałów
f1 = 1001.2; fs1 = 8000;
f2 = 303.1;  fs2 = 32000;
f3 = 2110.4; fs3 = 48000;

%% 2. Generacja sygnałów sinusoidalnych
t1 = (0:1/fs1:t-1/fs1)'; x1 = sin(2*pi*f1*t1);
t2 = (0:1/fs2:t-1/fs2)'; x2 = sin(2*pi*f2*t2);
t3 = (0:1/fs3:t-1/fs3)'; x3 = sin(2*pi*f3*t3);
t_target = (0:1/fs_target:t-1/fs_target)';
x4_expected = sin(2*pi*f1*t_target) + sin(2*pi*f2*t_target) + sin(2*pi*f3*t_target);

%% 3. Repróbkowanie metodą resample
[p1, q1] = rat(fs_target/fs1);
[p2, q2] = rat(fs_target/fs2);
x1_resampled = resample(x1, p1, q1);
x2_resampled = resample(x2, p2, q2);
x3_resampled = x3;

min_len = min([length(x1_resampled), length(x2_resampled), length(x3_resampled)]);
x1_resampled = x1_resampled(1:min_len);
x2_resampled = x2_resampled(1:min_len);
x3_resampled = x3_resampled(1:min_len);
x4_resampled = x1_resampled + x2_resampled + x3_resampled;
x4_expected_resampled = x4_expected(1:min_len);
mse_resample = mean((x4_resampled - x4_expected_resampled).^2);

%% 4. Repróbkowanie metodą ręczną (upsampling + filtracja + decymacja)

% Sygnał x1
[p1, q1] = rat(fs_target/fs1);
x1_up = upsample(x1, p1);  % Nadpróbkowanie (wstawianie zer)
cutoff1 = 1 / max(p1, q1);  % Częstotliwość odcięcia filtru (w normalizowanych jednostkach)
h1 = fir1(256, cutoff1);   % Filtr interpolacyjny i antyaliasingowy
delay1 = floor(length(h1)/2);
x1_filt = conv(x1_up, h1, 'full');       % Filtracja
x1_filt = x1_filt(delay1+1:end-delay1);  % Kompensacja opóźnienia
x1_manual = downsample(x1_filt, q1);     % Decymacja (wybranie co q1-tej próbki)

% Sygnał x2
[p2, q2] = rat(fs_target/fs2);
x2_up = upsample(x2, p2);
cutoff2 = 1 / max(p2, q2);
h2 = fir1(256, cutoff2);
delay2 = floor(length(h2)/2);
x2_filt = conv(x2_up, h2, 'full');
x2_filt = x2_filt(delay2+1:end-delay2);
x2_manual = downsample(x2_filt, q2);

% Sygnał x3 (już w fs_target, więc bez zmian)
x3_manual = x3;

% Dopasowanie długości
min_len_manual = min([length(x1_manual), length(x2_manual), length(x3_manual), length(x4_expected)]);
x1_manual = x1_manual(1:min_len_manual);
x2_manual = x2_manual(1:min_len_manual);
x3_manual = x3_manual(1:min_len_manual);
x4_manual = x1_manual + x2_manual + x3_manual;
x4_expected_manual = x4_expected(1:min_len_manual);

% Obliczenie błędu MSE
mse_manual = mean((x4_manual - x4_expected_manual).^2);


%% 5. Repróbkowanie metodą analityczną (idealna wiedza o funkcji)
x1_analytical = sin(2*pi*f1*t_target);
x2_analytical = sin(2*pi*f2*t_target);
x3_analytical = sin(2*pi*f3*t_target);
x4_analytical = x1_analytical + x2_analytical + x3_analytical;

x4_expected_analytical = x4_expected; % pełna długość
mse_analytical = mean((x4_analytical - x4_expected_analytical).^2);





%% 6. Wypisanie błędów
fprintf('MSE (resample): %.10e\n', mse_resample);
fprintf('MSE (manual):   %.10e\n', mse_manual);
fprintf('MSE (analytical): %.10e\n', mse_analytical);

%% 7. Wykresy porównawcze (sygnały)
samples_to_show = 500;
figure('Name', 'Porównanie metod repróbkowania', 'Position', [100, 100, 1000, 800]);

subplot(3,1,1);
plot(t_target(1:samples_to_show), x4_resampled(1:samples_to_show), 'b', ...
     t_target(1:samples_to_show), x4_expected_resampled(1:samples_to_show), 'k--');
title('Metoda resample'); legend('resample', 'oczekiwany'); grid on;

subplot(3,1,2);
plot(t_target(1:samples_to_show), x4_manual(1:samples_to_show), 'g', ...
     t_target(1:samples_to_show), x4_expected_manual(1:samples_to_show), 'k--');
title('Metoda ręczna (upsample + filtr + downsample)'); legend('manual', 'oczekiwany'); grid on;

subplot(3,1,3);
plot(t_target(1:samples_to_show), x4_analytical(1:samples_to_show), 'r', ...
     t_target(1:samples_to_show), x4_expected_analytical(1:samples_to_show), 'k--');
title('Metoda analityczna'); legend('analityczny', 'oczekiwany'); grid on;

%% 8. Analiza widmowa
NFFT = 2^nextpow2(length(x4_expected));
f = fs_target/2 * linspace(0, 1, NFFT/2+1);

X_resample = fft(x4_resampled, NFFT) / length(x4_resampled);
X_manual = fft(x4_manual, NFFT) / length(x4_manual);
X_analytical = fft(x4_analytical, NFFT) / length(x4_analytical);

figure('Name', 'Widma sygnałów dla różnych metod', 'Position', [100, 100, 1000, 500]);
plot(f, 2*abs(X_resample(1:NFFT/2+1)), 'b');
hold on;
plot(f, 2*abs(X_manual(1:NFFT/2+1)), 'g');
plot(f, 2*abs(X_analytical(1:NFFT/2+1)), 'r');
title('Widmo sygnału wyjściowego'); xlabel('Częstotliwość [Hz]'); ylabel('Amplituda');
legend('resample', 'manual', 'analytical'); grid on;
xlim([0, 3000]);

%% 8. Miksowanie plików WAV
try
    [wav1, fs_wav1] = audioread('x1.wav');
    [wav2, fs_wav2] = audioread('x2.wav');
    if size(wav1,2) > 1, wav1 = mean(wav1,2); end
    if size(wav2,2) > 1, wav2 = mean(wav2,2); end

    [p1, q1] = rat(fs_target/fs_wav1);
    [p2, q2] = rat(fs_target/fs_wav2);
    wav1 = resample(wav1, p1, q1);
    wav2 = resample(wav2, p2, q2);

    min_wav_len = min(length(wav1), length(wav2));
    wav_mixed = wav1(1:min_wav_len) + wav2(1:min_wav_len);
    wav_mixed = wav_mixed / max(abs(wav_mixed));

    audiowrite('mixed_48kHz.wav', wav_mixed, fs_target);
    sound(wav_mixed, fs_target);
    pause();
catch
    warning('Nie znaleziono plików x1.wav lub x2.wav.');
end
