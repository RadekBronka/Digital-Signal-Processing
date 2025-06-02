% Wczytanie i przygotowanie
[x, fs] = audioread('DontWorryBeHappy.wav');
x = mean(x, 2); % mono
N = 128; K = N/2;
x = x(1:floor(length(x)/K)*K); % dopasuj długość
w = sin(pi * ((0:N-1)+0.5) / N)';

% Macierz MDCT
A = sqrt(4/N) * cos(2*pi/N * ((0:K-1)' + 0.5) * ((0:N-1) + 0.5 + N/4));

% Analiza
numFrames = length(x)/K - 1;
Xmdct = zeros(K, numFrames);
for i = 1:numFrames
    idx = (i-1)*K + 1 : (i-1)*K + N;
    frame = x(idx) .* w;
    Xmdct(:, i) = A * frame;
end
Xmax = max(abs(Xmdct(:)));

% Szukaj Q dla 64 kbps
target_bitrate = 64000; % bps
Q_target = ceil(target_bitrate / (fs/2));

fprintf("Minimalna liczba bitów Q dla 64 kbps: %d\n", Q_target);

% Przetestuj różne Q
Q_values = Q_target:16; % od Q_target do 16 bitów
mse_values = zeros(size(Q_values));

for j = 1:length(Q_values)
    Q = Q_values(j);
    levels = 2^Q;
    Xq = round(Xmdct / Xmax * (levels/2 - 1));
    Xrec = Xq * Xmax / (levels/2 - 1);

    % Synteza
    x_rec = zeros(K * (numFrames + 1), 1);
    for i = 1:numFrames
        frame_rec = A' * Xrec(:, i);
        idx = (i-1)*K + 1 : (i-1)*K + N;
        x_rec(idx) = x_rec(idx) + frame_rec .* w;
    end
    x_rec = x_rec(1:length(x));
    mse_values(j) = mean((x - x_rec).^2);
end

% Wykres
figure;
plot(Q_values, 10*log10(mse_values), '-o');
xlabel('Liczba bitów Q');
ylabel('MSE [dB]');
title('Jakość rekonstrukcji vs liczba bitów Q');
grid on;
