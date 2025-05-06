% --- Parametry ---
fs1 = 8000;    % fs dla x1
fs2 = 32000;   % fs dla x2
fs3 = 48000;   % docelowe fs
t1 = 0:1/fs1:1-1/fs1;
t2 = 0:1/fs2:1-1/fs2;
t3 = 0:1/fs3:1-1/fs3; % już w docelowym

f1 = 1001.2;
f2 = 303.1;
f3 = 2110.4;

% --- Generowanie sygnałów ---
x1 = sin(2*pi*f1*t1);
x2 = sin(2*pi*f2*t2);
x3 = sin(2*pi*f3*t3); % już w fs3

%% --- Repróbkowanie x1: 8000 -> 48000 (L = 6, M = 1) ---
L1 = 6; M1 = 1;

% Upsampling przez L1 (dodanie zer)
x1_up = zeros(1, length(x1)*L1);
x1_up(1:L1:end) = x1;

% Interpolacja: filtr dolnoprzepustowy (np. FIR)
Nfir = 64;  % długość filtru
h1 = fir1(Nfir, 1/L1);  % filtr z pasmem 1/L
x1_filt = filter(h1, 1, x1_up);  % filtracja (interpolacja)

% Decymacja (tu niepotrzebna, bo M=1)
x1_resampled = x1_filt;

%% --- Repróbkowanie x2: 32000 -> 48000 (L = 3, M = 2) ---
L2 = 3; M2 = 2;

% Upsampling przez L2
x2_up = zeros(1, length(x2)*L2);
x2_up(1:L2:end) = x2;

% Interpolacja filtrem
h2 = fir1(Nfir, 1/max(L2, M2));
x2_filt = filter(h2, 1, x2_up);

% Decymacja przez M2
x2_resampled = x2_filt(1:M2:end);

%% --- x3 nie trzeba zmieniać ---
x3_resampled = x3;

%% --- Dopasowanie długości ---
min_len = min([length(x1_resampled), length(x2_resampled), length(x3_resampled)]);
x1_resampled = x1_resampled(1:min_len);
x2_resampled = x2_resampled(1:min_len);
x3_resampled = x3_resampled(1:min_len);

% --- Sumowanie sygnałów ---
x4 = x1_resampled + x2_resampled + x3_resampled;

% --- Analityczny wzorzec ---
t4 = 0:1/fs3:1-1/fs3;
x4_analytical = sin(2*pi*f1*t4) + sin(2*pi*f2*t4) + sin(2*pi*f3*t4);
x4_analytical = x4_analytical(1:min_len);

%% --- Porównanie ---
figure;
plot(t4(1:min_len), x4, 'b', t4(1:min_len), x4_analytical, 'r--');
legend('x4 (manualne resampling)', 'x4 (analityczny)');
xlabel('Czas [s]');
ylabel('Amplituda');
title('Porównanie ręcznie złożonego sygnału z analitycznym');

%% --- Odsłuch ---
disp('Odtwarzanie x4 (manualne przetwarzanie)...');
sound(x4, fs3);
pause();
disp('Odtwarzanie x4_analytical...');
sound(x4_analytical, fs3);
pause();


%% --- Miksowanie rzeczywistych plików WAV (x1.wav i x2.wav) ---
% Wczytaj dwa pliki
[x1_real, fs1_real] = audioread('x1.wav');
[x2_real, fs2_real] = audioread('x2.wav');

% Konwersja stereo -> mono jeśli trzeba
if size(x1_real, 2) > 1
    x1_real = mean(x1_real, 2);
end
if size(x2_real, 2) > 1
    x2_real = mean(x2_real, 2);
end

% Repróbkowanie do fs3
x1_real_resampled = resample(x1_real, fs3, fs1_real);
x2_real_resampled = resample(x2_real, fs3, fs2_real);

% Dopasowanie długości
min_len_real = min(length(x1_real_resampled), length(x2_real_resampled));
x1_real_resampled = x1_real_resampled(1:min_len_real);
x2_real_resampled = x2_real_resampled(1:min_len_real);

% Sumowanie
x_real_mix = x1_real_resampled + x2_real_resampled;

% Normalizacja (by uniknąć przesterowania)
x_real_mix = x_real_mix / max(abs(x_real_mix));

% Odsłuch
disp('Odtwarzanie zmiksowanego sygnału z plików WAV...');
sound(x_real_mix, fs3);

% Zapis do pliku
audiowrite('mix_real.wav', x_real_mix, fs3);

% Wykres
t_real = (0:min_len_real-1)/fs3;
figure;
plot(t_real, x_real_mix);
xlabel('Czas [s]');
ylabel('Amplituda');
title('Zmiksowany sygnał z plików x1.wav i x2.wav');
