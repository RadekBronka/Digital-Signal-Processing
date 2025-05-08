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

% Dopasowanie długości
min_len = min([length(x1_resampled), length(x2_resampled), length(x3_resampled)]);
x1_resampled = x1_resampled(1:min_len);
x2_resampled = x2_resampled(1:min_len);
x3_resampled = x3_resampled(1:min_len);
x4_resampled = x1_resampled + x2_resampled + x3_resampled;
x4_expected = x4_expected(1:min_len);

%% 4. Obliczenie błędu MSE
mse = mean((x4_resampled - x4_expected).^2);
fprintf('MSE między sygnałem wynikowym a oczekiwanym: %.10e\n', mse);

%% 5. Wykresy porównawcze
samples_to_show = 500;

figure('Name', 'Porównanie repróbkowanych sygnałów', 'Position', [100, 100, 900, 700]);
subplot(3,1,1);
plot(t_target(1:samples_to_show), x1_resampled(1:samples_to_show), 'r');
title('Sygnał x1 po resample'); grid on;

subplot(3,1,2);
plot(t_target(1:samples_to_show), x2_resampled(1:samples_to_show), 'g');
title('Sygnał x2 po resample'); grid on;

subplot(3,1,3);
plot(t_target(1:samples_to_show), x4_resampled(1:samples_to_show), 'b', ...
     t_target(1:samples_to_show), x4_expected(1:samples_to_show), 'k--');
title('Sygnał wynikowy vs. oczekiwany'); legend('Wynikowy', 'Oczekiwany'); grid on;

%% 6. Analiza widmowa
figure('Name', 'Widmo sygnału wynikowego', 'Position', [100, 100, 900, 500]);
NFFT = 2^nextpow2(length(x4_resampled));
X4 = fft(x4_resampled, NFFT) / length(x4_resampled);
f = fs_target/2 * linspace(0, 1, NFFT/2+1);
plot(f, 2*abs(X4(1:NFFT/2+1)));
title('Widmo sygnału x4'); xlabel('Hz'); ylabel('Amplituda'); grid on; xlim([0, 3000]);
hold on;
line([f1 f1], [0 0.5], 'Color', 'r', 'LineStyle', '--');
line([f2 f2], [0 0.5], 'Color', 'g', 'LineStyle', '--');
line([f3 f3], [0 0.5], 'Color', 'b', 'LineStyle', '--');

%% 7. Odsłuch
x4_play = x4_resampled / max(abs(x4_resampled));
sound(x4_play, fs_target);
pause();

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
