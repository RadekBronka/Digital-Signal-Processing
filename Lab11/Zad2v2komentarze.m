% === Wczytanie pliku audio ===
[x_orig, fs] = audioread('DontWorryBeHappy.wav');  % Wczytaj plik audio i częstotliwość próbkowania
x_orig = mean(x_orig, 2);  % Zamień stereo na mono (uśrednienie kanałów)

% === Parametry testowe: różne długości transformacji MDCT ===
N_values = [32, 128];  % Długości ramki MDCT do porównania

% === Pętla po różnych długościach ramki N ===
for n_idx = 1:length(N_values)
    N = N_values(n_idx);  % Długość jednej ramki (przed transformacją MDCT)
    K = N / 2;            % Przesunięcie ramki — MDCT stosuje przesunięcie o połowę okna
    
    % === Dopasowanie sygnału do pełnej liczby ramek ===
    x = x_orig(1:floor(length(x_orig)/K)*K);  % Przycięcie sygnału do wielokrotności K
    
    % === Okno sinusoidalne (minimalizacja przecieków widmowych) ===
    n = 0:N-1;
    w = sin(pi * (n + 0.5) / N)';  % Okno o długości N; transpozycja, by był wektorem kolumnowym

    % === Liczba ramek MDCT ===
    numFrames = length(x)/K - 1;  % Ostatnia ramka niepełna — dlatego -1
    
    % === Macierz MDCT ===
    % Współczynniki kosinusowe — przygotowanie do transformacji MDCT
    A = sqrt(4/N) * cos(2*pi/N * ((0:K-1)' + 0.5) * ((0:N-1) + 0.5 + N/4));

    % === Transformacja MDCT (analiza) ===
    Xmdct = zeros(K, numFrames);  % Zainicjalizuj macierz na współczynniki MDCT
    for i = 1:numFrames
        idx = (i-1)*K + 1 : (i-1)*K + N;  % Indeksy aktualnej ramki
        frame = x(idx) .* w;             % Zastosuj okno do ramki
        Xmdct(:, i) = A * frame;         % Oblicz współczynniki MDCT dla ramki
    end

    % === Kwantyzacja 16-bitowa ===
    % Symulacja bezstratnej kwantyzacji — zakładamy bardzo wysoką precyzję
    Xmax = max(abs(Xmdct(:)));  % Maksymalna wartość amplitudy MDCT
    Xq = round(Xmdct / Xmax * (2^15 - 1));  % Kwantyzacja do zakresu [-32767, 32767]
    Xrec = Xq * Xmax / (2^15 - 1);          % Dekwantyzacja (przywrócenie wartości rzeczywistych)

    % === Rekonstrukcja (synteza) — odwrotna MDCT ===
    x_rec = zeros(K * (numFrames + 1), 1);  % Bufor na zrekonstruowany sygnał
    for i = 1:numFrames
        frame_rec = A' * Xrec(:, i);        % Odtwórz próbki ramki z MDCT
        idx = (i-1)*K + 1 : (i-1)*K + N;     % Indeksy miejsca, gdzie dodać ramkę
        x_rec(idx) = x_rec(idx) + frame_rec .* w;  % Dodaj ramkę z oknem (nakładanie się ramek)
    end

    % === Dopasowanie długości końcowej i normalizacja ===
    x_rec = x_rec(1:length(x));              % Upewnij się, że długość pasuje
    x_rec = x_rec / max(abs(x_rec));         % Normalizacja amplitudy do [-1, 1]

    % === Obliczenie błędu rekonstrukcji (MSE) ===
    mse = mean((x - x_rec).^2);              % Średni błąd kwadratowy (jakość rekonstrukcji)
    fprintf('N = %d | MSE rekonstrukcji (Q=16 bitów): %.2e\n', N, mse);  % Wypisz wynik

    % === Wykres porównujący fragmenty oryginalnego i zrekonstruowanego sygnału ===
    figure;
    t = (0:1000)/fs;  % Oś czasu dla 1000 próbek (~0.02s)
    plot(t, x(1:1001), 'b', t, x_rec(1:1001), 'r--');  % Porównanie fragmentu sygnału
    legend('Oryginalny', 'Zrekonstruowany');
    xlabel('Czas [s]'); ylabel('Amplituda');
    title(sprintf('N = %d | Porównanie sygnałów oryginalnego i zrekonstruowanego', N));

    % === Odtwarzanie dźwięku (odsłuch efektów kompresji i rekonstrukcji) ===
    fprintf('Odtwarzanie fragmentu sygnału zrekonstruowanego (N=%d)...\n', N);
    soundsc(x_rec, fs);  % Odtwórz zrekonstruowany sygnał z normalizacją głośności
    pause(length(x_rec)/fs + 1);  % Czekaj, aż sygnał się odtworzy (plus 1 sekunda przerwy)
end
