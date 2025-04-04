% Parametry
fp = 64e3;      % Hz
fs = 128e3;     % Hz
wp = 2*pi*fp;   % rad/s
ws = 2*pi*fs;   % rad/s
Rp = 3;         % maksymalne zafalowanie (dB)
Rs = 40;        % minimalne tłumienie (dB)

types = {'butter', 'cheby1', 'cheby2', 'ellip'};
titles = {'Butterworth', 'Czebyszew I', 'Czebyszew II', 'Eliptyczny'};
H_all = {};
orders = zeros(1,4);

for i = 1:4
    switch types{i}
        case 'butter'
            [n, wn] = buttord(wp, ws, Rp, Rs, 's');
            [b, a] = butter(n, wn, 's');
        case 'cheby1'
            [n, wn] = cheb1ord(wp, ws, Rp, Rs, 's');
            [b, a] = cheby1(n, Rp, wn, 's');
        case 'cheby2'
            [n, wn] = cheb2ord(wp, ws, Rp, Rs, 's');
            [b, a] = cheby2(n, Rs, wn, 's');
        case 'ellip'
            [n, wn] = ellipord(wp, ws, Rp, Rs, 's');
            [b, a] = ellip(n, Rp, Rs, wn, 's');
    end
    orders(i) = n;
    H_all{i} = tf(b, a);
end

% Wykres biegunów
figure('Name', 'Bieguny filtrów');
for i = 1:4
    subplot(2,2,i);
    pzmap(H_all{i});
    title(['Bieguny – ' titles{i} ', rząd = ' num2str(orders(i))]);
    grid on;
end

% Charakterystyki amplitudowe
f = linspace(10e3, 200e3, 1000);
w = 2*pi*f;

figure('Name', 'Charakterystyki amplitudowe');
hold on;
for i = 1:4
    [mag, ~] = bode(H_all{i}, w);
    mag = squeeze(mag);
    plot(f/1e3, 20*log10(mag), 'DisplayName', titles{i});
end
xlabel('Częstotliwość [kHz]');
ylabel('Wzmocnienie [dB]');
title('Porównanie charakterystyk filtrów');
legend show;
grid on;
