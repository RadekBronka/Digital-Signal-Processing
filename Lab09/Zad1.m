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
