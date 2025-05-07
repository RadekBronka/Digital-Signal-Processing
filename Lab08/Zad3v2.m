%% Zadanie 3 - Filtr interpolatora i decymatora, mikser audio
clear all;
close all;
clc;

%% 1. Definicja parametrów sygnałów
% Parametry sygnałów
fs_target = 48000; % Hz - docelowa częstotliwość próbkowania
t = 1; % s - czas trwania sygnału

% Sygnały składowe
f1 = 1001.2; % Hz
fs1 = 8000;  % Hz
f2 = 303.1;  % Hz
fs2 = 32000; % Hz
f3 = 2110.4; % Hz
fs3 = 48000; % Hz

%% 2. Generacja sygnałów sinusoidalnych
% Generacja wektorów czasu i sygnałów
t1 = (0:1/fs1:t-1/fs1)';
t2 = (0:1/fs2:t-1/fs2)';
t3 = (0:1/fs3:t-1/fs3)';
t_target = (0:1/fs_target:t-1/fs_target)';

x1 = sin(2*pi*f1*t1);
x2 = sin(2*pi*f2*t2);
x3 = sin(2*pi*f3*t3);

% Generacja teoretycznego sygnału wzorcowego
x4_expected = sin(2*pi*f1*t_target) + sin(2*pi*f2*t_target) + sin(2*pi*f3*t_target);

%% 3. Repróbkowanie sygnałów metodą 1 - użycie funkcji resample
% Metoda 1: Korzystamy z gotowej funkcji resample
% Funkcja resample wykonuje zarówno nadpróbkowanie, interpolację jak i decymację
disp('Metoda 1: Repróbkowanie za pomocą funkcji resample');

% Obliczenie współczynników interpolacji/decymacji
[p1, q1] = rat(fs_target/fs1); % p1=6, q1=1 (fs_target/fs1 = 48000/8000 = 6)
[p2, q2] = rat(fs_target/fs2); % p2=3, q2=2 (fs_target/fs2 = 48000/32000 = 3/2)

% Repróbkowanie sygnałów
x1_resampled = resample(x1, p1, q1);
x2_resampled = resample(x2, p2, q2); 
x3_resampled = x3; % już ma odpowiednią częstotliwość

%% 4. Repróbkowanie sygnałów metodą 2 - ręczna implementacja nadpróbkowania i interpolacji
disp('Metoda 2: Ręczna implementacja nadpróbkowania i interpolacji');

% Sygnał x1 (fs1 = 8000 Hz → fs_target = 48000 Hz)
% Współczynnik nadpróbkowania = 6
up_factor1 = fs_target / fs1;

% Nadpróbkowanie przez wstawienie zer (upsampling)
x1_upsampled = zeros(length(x1) * up_factor1, 1);
x1_upsampled(1:up_factor1:end) = x1;

% Projektowanie filtru dolnoprzepustowego jako filtru interpolującego
% Częstotliwość odcięcia = fs1/2 / fs_target = 4000/48000 = 1/12
cutoff1 = (fs1/2) / fs_target;
filter_order1 = 48; % Rząd filtru
b_interp1 = fir1(filter_order1, cutoff1, 'low') * up_factor1;

% Filtracja interpolująca
x1_manual = filter(b_interp1, 1, x1_upsampled);

% Kompensacja opóźnienia filtru
delay1 = filter_order1 / 2;
x1_manual = x1_manual(delay1+1:end);
x1_manual = [x1_manual; zeros(delay1, 1)];

% Sygnał x2 (fs2 = 32000 Hz → fs_target = 48000 Hz)
% Współczynnik nadpróbkowania = 3/2 = 3
% Współczynnik decymacji = 2
up_factor2 = 3;
down_factor2 = 2;

% Nadpróbkowanie przez wstawienie zer
x2_upsampled = zeros(length(x2) * up_factor2, 1);
x2_upsampled(1:up_factor2:end) = x2;

% Filtr interpolujący
cutoff2 = (fs2/2) / (fs2*up_factor2);
filter_order2 = 48;
b_interp2 = fir1(filter_order2, cutoff2, 'low') * up_factor2;

% Filtracja interpolująca
x2_filtered = filter(b_interp2, 1, x2_upsampled);

% Kompensacja opóźnienia
delay2 = filter_order2 / 2;
x2_filtered = x2_filtered(delay2+1:end);
x2_filtered = [x2_filtered; zeros(delay2, 1)];

% Decymacja (wybranie co down_factor2-tej próbki)
x2_manual = x2_filtered(1:down_factor2:end);

% Sygnał x3 już ma docelową częstotliwość próbkowania
x3_manual = x3;

%% 5. Miksowanie sygnałów - dopasowanie długości i dodawanie

% Znajdź minimalną długość dla obu metod
min_length1 = min([length(x1_resampled), length(x2_resampled), length(x3_resampled)]);
min_length2 = min([length(x1_manual), length(x2_manual), length(x3_manual)]);

% Przytnij sygnały do jednakowej długości
% Metoda 1 - resample
x1_resampled = x1_resampled(1:min_length1);
x2_resampled = x2_resampled(1:min_length1);
x3_resampled = x3_resampled(1:min_length1);
x4_resampled = x1_resampled + x2_resampled + x3_resampled;

% Metoda 2 - ręczna implementacja
x1_manual = x1_manual(1:min_length2);
x2_manual = x2_manual(1:min_length2);
x3_manual = x3_manual(1:min_length2);
x4_manual = x1_manual + x2_manual + x3_manual;

% Przytnij sygnał oczekiwany do porównania
x4_expected_1 = x4_expected(1:min_length1);
x4_expected_2 = x4_expected(1:min_length2);

%% 6. Obliczenie błędów MSE dla podstawowych metod
mse_resampled = mean((x4_resampled - x4_expected_1).^2);
mse_manual = mean((x4_manual - x4_expected_2).^2);

% Wyświetlenie wyników
fprintf('Porównanie metod repróbkowania (MSE):\n');
fprintf('1. Funkcja resample: %.10e\n', mse_resampled);
fprintf('2. Ręczna implementacja: %.10e\n', mse_manual);

%% 7. Wizualizacja wyników - porównanie repróbkowanych sygnałów

% Przygotowanie osi czasu do wyświetlania
samples_to_show = 500; % Liczba próbek do wyświetlenia

figure('Name', 'Porównanie metod repróbkowania', 'Position', [100, 100, 900, 700]);

% Porównanie sygnałów x1 po repróbkowaniu różnymi metodami
subplot(3,1,1);
plot(t_target(1:samples_to_show), x1_resampled(1:samples_to_show), 'r-', ...
     t_target(1:samples_to_show), x1_manual(1:samples_to_show), 'g--');
title('Porównanie repróbkowania sygnału x1');
legend('resample', 'Ręczna implementacja');
grid on;

% Porównanie sygnałów x2 po repróbkowaniu różnymi metodami
subplot(3,1,2);
plot(t_target(1:samples_to_show), x2_resampled(1:samples_to_show), 'r-', ...
     t_target(1:samples_to_show), x2_manual(1:samples_to_show), 'g--');
title('Porównanie repróbkowania sygnału x2');
legend('resample', 'Ręczna implementacja');
grid on;

% Porównanie sygnałów wynikowych z analitycznym x4
subplot(3,1,3);
plot(t_target(1:samples_to_show), x4_resampled(1:samples_to_show), 'r-', ...
     t_target(1:samples_to_show), x4_manual(1:samples_to_show), 'g--', ...
     t_target(1:samples_to_show), x4_expected_1(1:samples_to_show), 'k--');
title('Porównanie sygnału wynikowego x4 z różnych metod');
legend('resample', 'Ręczna implementacja', 'Analityczny');
grid on;

%% 8. Analiza widmowa dla najlepszej metody (resample)
figure('Name', 'Analiza widmowa', 'Position', [100, 100, 900, 500]);

% Widmo sygnału wynikowego z metody resample
NFFT = 2^nextpow2(length(x4_resampled));
X4 = fft(x4_resampled, NFFT) / length(x4_resampled);
f = fs_target/2 * linspace(0, 1, NFFT/2+1);
plot(f, 2*abs(X4(1:NFFT/2+1)));
title('Widmo sygnału wynikowego x4 (metoda resample)');
xlabel('Częstotliwość [Hz]');
ylabel('Amplituda');
grid on;
xlim([0, 3000]); % Ograniczenie zakresu dla lepszej widoczności
% Zaznaczenie częstotliwości składowych sygnałów
hold on;
line([f1 f1], [0 0.5], 'Color', 'r', 'LineStyle', '--');
line([f2 f2], [0 0.5], 'Color', 'g', 'LineStyle', '--');
line([f3 f3], [0 0.5], 'Color', 'b', 'LineStyle', '--');
legend('Widmo sygnału', ['f1 = ' num2str(f1) ' Hz'], ['f2 = ' num2str(f2) ' Hz'], ['f3 = ' num2str(f3) ' Hz']);

%% 9. Odsłuch sygnałów
% Normalizacja sygnałów przed odsłuchem
x1_norm = x1_resampled/max(abs(x1_resampled));
x2_norm = x2_resampled/max(abs(x2_resampled));
x3_norm = x3_resampled/max(abs(x3_resampled));
x4_norm = x4_resampled/max(abs(x4_resampled));

% Odsłuch sygnałów składowych
disp('Odsłuch sygnału x1:');
sound(x1_norm, fs_target);
pause(1.5); % Krótka pauza między sygnałami

disp('Odsłuch sygnału x2:');
sound(x2_norm, fs_target);
pause(1.5);

disp('Odsłuch sygnału x3:');
sound(x3_norm, fs_target);
pause(1.5);

% Odsłuch sygnału wynikowego
disp('Odsłuch sygnału wynikowego (resample):');
sound(x4_norm, fs_target);

%% 10. Wczytywanie i miksowanie plików WAV
% Wczytywanie plików WAV
try
    [wav1, fs_wav1] = audioread('x1.wav');
    [wav2, fs_wav2] = audioread('x2.wav');
    
    % Konwersja do mono jeśli stereo
    if size(wav1, 2) > 1
        wav1 = mean(wav1, 2);
    end
    if size(wav2, 2) > 1
        wav2 = mean(wav2, 2);
    end
    
    % Repróbkowanie do 48000 Hz (metoda resample, jako najlepsza)
    [p_wav1, q_wav1] = rat(fs_target/fs_wav1); 
    [p_wav2, q_wav2] = rat(fs_target/fs_wav2);
    
    wav1_resampled = resample(wav1, p_wav1, q_wav1);
    wav2_resampled = resample(wav2, p_wav2, q_wav2);
    
    % Dostosowanie długości sygnałów
    min_length_wav = min(length(wav1_resampled), length(wav2_resampled));
    wav1_resampled = wav1_resampled(1:min_length_wav);
    wav2_resampled = wav2_resampled(1:min_length_wav);
    
    % Miksowanie sygnałów
    wav_mixed = wav1_resampled + wav2_resampled;
    
    % Normalizacja
    wav_mixed = wav_mixed / max(abs(wav_mixed));
    
    % Zapis do pliku
    audiowrite('mixed_48kHz.wav', wav_mixed, fs_target);
    
    % Odsłuch
    disp('Odsłuch zmiksowanych plików WAV (48kHz):');
    sound(wav_mixed, fs_target);
    
catch
    warning('Nie udało się wczytać plików WAV. Upewnij się, że pliki x1.wav i x2.wav znajdują się w katalogu roboczym.');
end
