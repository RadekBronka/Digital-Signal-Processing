% Wczytanie pliku
[x_orig, fs] = audioread('DontWorryBeHappy.wav');
x_orig = mean(x_orig, 2);  % Mono

% Testowane wartości N
N_values = [32, 128];

for n_idx = 1:length(N_values)
    N = N_values(n_idx);
    K = N / 2;
    
    % Dopasowanie długości do pełnych ramek
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
    fprintf('N = %d | MSE rekonstrukcji (Q=16 bitów): %.2e\n', N, mse);
    
     % Wykres porównawczy fragmentu
    figure;
    t = (0:1000)/fs;
    plot(t, x(1:1001), 'b', t, x_rec(1:1001), 'r--');
    legend('Oryginalny', 'Zrekonstruowany');
    xlabel('Czas [s]'); ylabel('Amplituda');
    title(sprintf('N = %d | Porównanie sygnałów oryginalnego i zrekonstruowanego', N));
   
    
    % Odsłuch fragmentu zrekonstruowanego
    fprintf('Odtwarzanie fragmentu sygnału zrekonstruowanego (N=%d)...\n', N);
    soundsc(x_rec(), fs);
    pause(length(x_rec())/fs + 1);
end
