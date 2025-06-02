% Wczytanie pliku
[x, fs] = audioread('DontWorryBeHappy.wav');
x = mean(x, 2);              % Mono
N = 128;                     % Długość ramki MDCT
K = N/2;                     % Liczba współczynników

% Dopasowanie długości sygnału (całkowita liczba ramek)
x = x(1:floor(length(x)/K)*K);

% Okno sinusoidalne
n = 0:N-1;
w = sin(pi * (n + 0.5) / N)';

% Liczba ramek
numFrames = length(x)/K - 1;

% Analiza MDCT
A = sqrt(4/N) * cos(2*pi/N * ((0:K-1)' + 0.5) * ((0:N-1) + 0.5 + N/4));
Xmdct = zeros(K, numFrames);
for i = 1:numFrames
    idx = (i-1)*K + 1 : (i-1)*K + N;
    frame = x(idx) .* w;
    Xmdct(:, i) = A * frame;
end

% Kwantyzacja 16-bitowa (symulacja bezstratności)
Xmax = max(abs(Xmdct(:)));
Xq = round(Xmdct / Xmax * (2^15 - 1));
Xrec = Xq * Xmax / (2^15 - 1);  % Dekwantyzacja

% Synteza MDCT
x_rec = zeros(K * (numFrames + 1), 1);
for i = 1:numFrames
    frame_rec = A' * Xrec(:, i);
    idx = (i-1)*K + 1 : (i-1)*K + N;
    x_rec(idx) = x_rec(idx) + frame_rec .* w;
end

% Dopasowanie długości i normalizacja
x_rec = x_rec(1:length(x));
x_rec = x_rec / max(abs(x_rec));

% Obliczenie MSE
mse = mean((x - x_rec).^2);
fprintf('MSE rekonstrukcji (Q=16 bitów): %.2e\n', mse);

% Odsłuch (opcjonalnie)
% sound(x_rec, fs);

% Wykres porównawczy (fragment)
t = (0:1000)/fs;
figure;
plot(t, x(1:1001), 'b', t, x_rec(1:1001), 'r--');
legend('Oryginalny', 'Zrekonstruowany');
xlabel('Czas [s]'); ylabel('Amplituda');
title('Porównanie sygnałów oryginalnego i zrekonstruowanego');
