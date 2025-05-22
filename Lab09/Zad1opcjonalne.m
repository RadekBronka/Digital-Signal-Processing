clear;
close all;

% Parametry
fs = 8000;
t = 0:1/fs:1-1/fs;
SNRs = [10, 20, 40];
M  = 48;
mi = 0.008;

%% 1. SYGNAŁ HARMONICZNY
A1 = -0.5; f1 = 34.2;
A2 = 1.0;  f2 = 115.5;
dref = A1 * sin(2*pi*f1*t) + A2 * sin(2*pi*f2*t);

figure;
for i = 1:length(SNRs)
    SNR = SNRs(i);
    d = awgn(dref, SNR, 'measured');
    x = [d(1), d(1:end-1)];

    y = zeros(size(d));
    e = zeros(size(d));
    bx = zeros(M,1);
    h = zeros(M,1);

    for n = 1:length(x)
        bx = [x(n); bx(1:M-1)];
        y(n) = h' * bx;
        e(n) = d(n) - y(n);
        h = h + mi * e(n) * bx;
    end

    % Oblicz SNR po odszumieniu
    sygnal = mean(dref.^2);
    szum = mean((dref - y).^2);
    SNR_post = 10 * log10(sygnal / szum);

    % Subplot dla SNR harmonicznego sygnału
    subplot(3, 3, i); 
    plot(t, dref, 'k', 'LineWidth', 1); 
    title(['Syg. Oryg. SNR = ', num2str(SNR), ' dB']);
    xlabel('Czas [s]');
    ylabel('Amplituda');
    grid on;
    legend('Oryginalny');

    subplot(3, 3, i+3); 
    plot(t, d, 'r'); 
    title(['Syg. Zaszumiony SNR = ', num2str(SNR), ' dB']);
    xlabel('Czas [s]');
    ylabel('Amplituda');
    grid on;
    legend('Zaszumiony');

    subplot(3, 3, i+6); 
    plot(t, y, 'b'); 
    title(['Syg. Odszumiony SNR Po = ', num2str(SNR_post, '%.2f'), ' dB']);
    xlabel('Czas [s]');
    ylabel('Amplituda');
    grid on;
    legend('Odszumiony');
end

%% 2. SFM (ograniczone do 50 ms)
fc = 1000;
df = 500;
fm = 0.25;
dref = sin(2*pi*fc*t + (df/fm)*sin(2*pi*fm*t));

Nvis = round(0.05 * fs);  % 50 ms = 400 próbek
tt = t(1:Nvis);

figure;
for i = 1:length(SNRs)
    SNR = SNRs(i);
    d = awgn(dref, SNR, 'measured');
    x = [d(1), d(1:end-1)];

    y = zeros(size(d));
    e = zeros(size(d));
    bx = zeros(M,1);
    h = zeros(M,1);

    for n = 1:length(x)
        bx = [x(n); bx(1:M-1)];
        y(n) = h' * bx;
        e(n) = d(n) - y(n);
        h = h + mi * e(n) * bx;
    end

    % Oblicz SNR po odszumieniu
    sygnal = mean(dref(1:Nvis).^2);
    szum = mean((dref(1:Nvis) - y(1:Nvis)).^2);
    SNR_post = 10 * log10(sygnal / szum);

    % Subplot dla SFM
    subplot(3, 3, i); 
    plot(tt, dref(1:Nvis), 'k', 'LineWidth', 1); 
    title(['Syg. Oryg. SNR = ', num2str(SNR), ' dB']);
    xlabel('Czas [s]');
    ylabel('Amplituda');
    grid on;
    legend('Oryginalny');

    subplot(3, 3, i+3); 
    plot(tt, d(1:Nvis), 'r'); 
    title(['Syg. Zaszumiony SNR = ', num2str(SNR), ' dB']);
    xlabel('Czas [s]');
    ylabel('Amplituda');
    grid on;
    legend('Zaszumiony');

    subplot(3, 3, i+6); 
    plot(tt, y(1:Nvis), 'b'); 
    title(['Syg. Odszumiony SNR Po = ', num2str(SNR_post, '%.2f'), ' dB']);
    xlabel('Czas [s]');
    ylabel('Amplituda');
    grid on;
    legend('Odszumiony');
end


%% 3. NAGRANIE AUDIO
[audio, fs] = audioread('alfaromeo164.wav');
if size(audio,2) > 1
    audio = audio(:,1); 
end
audio = audio(1:min(end, fs)); % Przytnij do 1 sekundy
t = (0:length(audio)-1)/fs;
dref = audio(:)';

figure;
for i = 1:length(SNRs)
    SNR = SNRs(i);
    d = awgn(dref, SNR, 'measured');
    x = [d(1), d(1:end-1)];

    y = zeros(size(d));
    e = zeros(size(d));
    bx = zeros(M,1);
    h = zeros(M,1);

    for n = 1:length(x)
        bx = [x(n); bx(1:M-1)];
        y(n) = h' * bx;
        e(n) = d(n) - y(n);
        h = h + mi * e(n) * bx;
    end

    % Oblicz SNR po odszumieniu
    sygnal = mean(dref.^2);
    szum = mean((dref - y).^2);
    SNR_post = 10 * log10(sygnal / szum);

    % Subplot dla audio
    subplot(3, 3, i); 
    plot(t, dref(1:min(end, fs)), 'k', 'LineWidth', 1); 
    title(['Syg. Oryg. SNR = ', num2str(SNR), ' dB']);
    xlabel('Czas [s]');
    ylabel('Amplituda');
    grid on;
    legend('Oryginalny');

    subplot(3, 3, i+3); 
    plot(t, d(1:min(end, fs)), 'r'); 
    title(['Syg. Zaszumiony SNR = ', num2str(SNR), ' dB']);
    xlabel('Czas [s]');
    ylabel('Amplituda');
    grid on;
    legend('Zaszumiony');

    subplot(3, 3, i+6); 
    plot(t, y(1:min(end, fs)), 'b'); 
    title(['Syg. Odszumiony SNR Po = ', num2str(SNR_post, '%.2f'), ' dB']);
    xlabel('Czas [s]');
    ylabel('Amplituda');
    grid on;
    legend('Odszumiony');
end

%% 4. OPCJONALNIE: SYGNAŁ MOWY
% Wczytanie
[sprawc, fs_sp] = audioread('mowa8000.wav');
if size(sprawc,2) > 1
    sprawc = sprawc(:,1);
end
% Upewnij się, że fs_sp == 8000
if fs_sp ~= 8000
    error('Plik mowa8000.wav musi mieć fs = 8000 Hz');
end

% Przytnij do 1 sekundy
dref = sprawc()';

% Parametry zaszumienia i adaptacji
SNR_in = 20;                  % wybierz wygodną wartość
M_list = [32, 48, 64, 96];    % długości filtru do przetestowania
mu_list = [0.005, 0.008, 0.015, 0.02];  % kroki adaptacji

best_SNR = -Inf;
best_M = NaN;
best_mu = NaN;
best_y = [];

% Grid search
for M = M_list
  for mu = mu_list
    d = awgn(dref, SNR_in, 'measured');
    x = [d(1), d(1:end-1)];
    
    y = zeros(size(d));
    h = zeros(M,1);
    bx = zeros(M,1);
    
    for n = 1:length(d)
      bx = [x(n); bx(1:end-1)];
      y(n) = h' * bx;
      e = d(n) - y(n);
      h = h + mu * e * bx;
    end
    
    % SNR po odszumieniu
    SNR_out = 10*log10( mean(dref.^2) / mean((dref - y).^2) );
    if SNR_out > best_SNR
      best_SNR = SNR_out;
      best_M = M;
      best_mu = mu;
      best_y = y;
      best_h = h;
    end
  end
end

fprintf('Najlepsze parametry: M = %d, mu = %.4f → SNR_out = %.2f dB\n', ...
        best_M, best_mu, best_SNR);

% Odsłuch do oceny jakości
disp('Odsłuchaj oryginał, zaszumiony i odszumiony:');
soundsc(dref, fs_sp);    pause();
soundsc(awgn(dref,SNR_in,'measured'), fs_sp);  pause();
soundsc(best_y, fs_sp);

% Wybieramy jakiś ostatni fragment (np. 100 ms) do analizy widmowej
Nfrag = round(fs_sp);
frag_ref   = dref(end-Nfrag+1:end);
frag_deno  = best_y(end-Nfrag+1:end);

% Charakterystyka filtru adaptacyjnego
figure;
freqz(best_h,1,1024,fs_sp);
title('Charakterystyka amplitudowo-częstotliwościowa filtru h');

% Widmo gęstości mocy końcowego fragmentu odszumionego
figure;
pwelch(frag_deno, hamming(256), 128, 1024, fs_sp);
title('Widmo gęstości mocy ostatniego fragmentu odszumionego sygnału');
