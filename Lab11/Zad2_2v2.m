% --- Parametry ---
[x_orig, fs] = audioread('DontWorryBeHappy.wav');
x_orig = mean(x_orig, 2); % Mono

N = 128; % Długość ramki MDCT
K = N/2;

% Dopasowanie długości sygnału do pełnych ramek
x = x_orig(1:floor(length(x_orig)/K)*K);

% Okno sinusoidalne
n = 0:N-1;
w = sin(pi * (n + 0.5) / N)';

% Liczba ramek
numFrames = length(x)/K - 1;

% Macierz MDCT
A = sqrt(4/N) * cos(2*pi/N * ((0:K-1)' + 0.5) * ((0:N-1) + 0.5 + N/4));

% Analiza MDCT
Xmdct = zeros(K, numFrames);
for i = 1:numFrames
    idx = (i-1)*K + 1 : (i-1)*K + N;
    frame = x(idx) .* w;
    Xmdct(:, i) = A * frame;
end

% Maksymalna wartość współczynników do normalizacji
Xmax = max(abs(Xmdct(:)));

% Lista Q do testów (liczba bitów kwantyzacji)
Q_list = 1:32;

mse_list = zeros(size(Q_list));
bitrate_list = zeros(size(Q_list));

for q_idx = 1:length(Q_list)
    Q = Q_list(q_idx);
    
    % Kwantyzacja Q-bitowa (symulacja)
    max_level = 2^(Q-1) - 1;
    Xq = round(Xmdct / Xmax * max_level);
    Xrec = Xq * Xmax / max_level; % Dekwantyzacja
    
    % Synteza MDCT
    x_rec = zeros(K*(numFrames+1), 1);
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
    mse_list(q_idx) = mse;
    
    % Obliczenie bitrate w bps i kbps
    bits_total = numFrames * K * Q;
    time_s = length(x) / fs;
    bitrate = bits_total / time_s;
    bitrate_list(q_idx) = bitrate / 1000; % kbps
    
    fprintf('Q = %2d bit | Bitrate = %7.2f kbps | MSE = %.3e\n', Q, bitrate_list(q_idx), mse);
end

% Wykres MSE vs Bitrate
figure;
plot(bitrate_list, mse_list, '-o', 'LineWidth', 2);
grid on;
xlabel('Bitrate [kbps]');
ylabel('MSE rekonstrukcji');
title('Jakość rekonstrukcji (MSE) vs przepływność (bitrate)');
set(gca,'YScale','log'); % lepsza wizualizacja różnic MSE
