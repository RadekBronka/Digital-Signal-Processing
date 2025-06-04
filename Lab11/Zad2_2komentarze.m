% --- Parametry ---
[x_orig, fs] = audioread('DontWorryBeHappy.wav');  % Wczytanie pliku audio
x_orig = mean(x_orig, 2);                          % Konwersja do mono przez uśrednienie kanałów

N = 128;           % Długość ramki MDCT
K = N/2;           % Przesunięcie ramki (typowe dla MDCT)

% --- Przycięcie sygnału do wielokrotności K ---
x = x_orig(1:floor(length(x_orig)/K)*K);  % Usunięcie nadmiarowych próbek z końca

% --- Okno sinusoidalne ---
n = 0:N-1;
w = sin(pi * (n + 0.5) / N)';  % Okno długości N, transponowane do wektora kolumnowego

% --- Liczba ramek (zachodzących na siebie) ---
numFrames = length(x)/K - 1;  % "-1", ponieważ MDCT potrzebuje 2K = N próbek na ramkę

% --- Macierz MDCT ---
% A to macierz przekształcenia dla transformacji MDCT
A = sqrt(4/N) * cos(2*pi/N * ((0:K-1)' + 0.5) * ((0:N-1) + 0.5 + N/4));

% --- Analiza MDCT ---
Xmdct = zeros(K, numFrames);  % Macierz na współczynniki MDCT
for i = 1:numFrames
    idx = (i-1)*K + 1 : (i-1)*K + N;     % Indeksy dla ramki
    frame = x(idx) .* w;                % Zastosowanie okna
    Xmdct(:, i) = A * frame;            % Transformacja MDCT
end

% --- Maksymalna wartość (do normalizacji przed kwantyzacją) ---
Xmax = max(abs(Xmdct(:)));  % Normalizacja pozwala efektywnie wykorzystać zakres bitowy

% --- Lista testowanych głębokości kwantyzacji (bitów na współczynnik) ---
Q_list = 1:32;

% --- Bufory wyników ---
mse_list = zeros(size(Q_list));         % Błędy rekonstrukcji (dla każdego Q)
bitrate_list = zeros(size(Q_list));     % Przepływność (dla każdego Q)

% === Pętla po wartościach Q ===
for q_idx = 1:length(Q_list)
    Q = Q_list(q_idx);  % Liczba bitów
    
    % --- Kwantyzacja współczynników MDCT ---
    max_level = 2^(Q-1) - 1;                        % Największa wartość w Q-bitowej skali symetrycznej
    Xq = round(Xmdct / Xmax * max_level);          % Kwantyzacja
    Xrec = Xq * Xmax / max_level;                  % Dekwantyzacja (przywrócenie wartości)

    % --- Synteza MDCT (odtwarzanie sygnału) ---
    x_rec = zeros(K*(numFrames+1), 1);             % Bufor na odtworzony sygnał
    for i = 1:numFrames
        frame_rec = A' * Xrec(:, i);               % Inwersja transformacji MDCT
        idx = (i-1)*K + 1 : (i-1)*K + N;            % Indeksy na które przypada ramka
        x_rec(idx) = x_rec(idx) + frame_rec .* w;  % Nakładanie ramek z oknem
    end

    % --- Dopasowanie długości i normalizacja ---
    x_rec = x_rec(1:length(x));                    % Ucięcie nadmiarowych próbek
    x_rec = x_rec / max(abs(x_rec));               % Normalizacja amplitudy [-1, 1]

    % --- Obliczenie MSE (średni błąd kwadratowy) ---
    mse = mean((x - x_rec).^2);                    % Miara jakości rekonstrukcji
    mse_list(q_idx) = mse;

    % --- Obliczenie bitrate (przepływności danych) ---
    bits_total = numFrames * K * Q;                % Liczba bitów zakodowanych współczynników
    time_s = length(x) / fs;                       % Czas trwania analizowanego sygnału
    bitrate = bits_total / time_s;                 % Bitrate [bps]
    bitrate_list(q_idx) = bitrate / 1000;          % Bitrate [kbps]

    % --- Informacja tekstowa ---
    fprintf('Q = %2d bit | Bitrate = %7.2f kbps | MSE = %.3e\n', Q, bitrate_list(q_idx), mse);
end

% === Wykres MSE vs bitrate ===
figure;
plot(bitrate_list, mse_list, '-o', 'LineWidth', 2);
grid on;
xlabel('Bitrate [kbps]');
ylabel('MSE rekonstrukcji');
title('Jakość rekonstrukcji (MSE) vs przepływność (bitrate)');
set(gca,'YScale','log');  % Skala logarytmiczna MSE — uwidacznia różnice jakości
