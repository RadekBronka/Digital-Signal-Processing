clc;
clear;

name = 'Radek';
f_c = 500;           % częstotliwość nośnej
f_pr = 16000;        % częstotliwość próbkowania
T = 0.1;             % czas trwania symbolu
t = 0:1/f_pr:T-1/f_pr;

% Zamiana znaków na bity
bits = reshape(dec2bin(double(name),8).',1,[]);

N_symbols = length(bits)/2;

signal = [];
for k = 1:N_symbols
    pair = bits(2*k-1:2*k);
    % Przypisanie fazy na podstawie pary bitów
    switch pair
        case '00'
            phase = 0;
        case '01'
            phase = pi/2;
        case '11'
            phase = pi;
        case '10'
            phase = 3*pi/2;
    end
    symbol = sin(2*pi*f_c*t + phase);
    signal = [signal symbol];
end

% Prezentacja fragmentu sygnału
czas_wyswietlania = 0.8;  
samples_to_plot = floor(czas_wyswietlania * f_pr);

figure;
plot((0:samples_to_plot-1)/f_pr, signal(1:samples_to_plot));
title(['QPSK sygnał dla imienia "' name '"']);
xlabel('Czas [s]');
ylabel('Amplituda');

% Demodulacja

recovered_bits = '';

% Tworzymy wzorce sinusoidalne dla czterech faz 
wave0 = sin(2*pi*f_c*t + 0);
wave90 = sin(2*pi*f_c*t + pi/2);
wave180 = sin(2*pi*f_c*t + pi);
wave270 = sin(2*pi*f_c*t + 3*pi/2);

for k = 1:N_symbols
    % długość jednego symbolu
    symbol_length = length(t);

    % od którego miejsca zaczyna się k-ty symbol
    start_index = (k-1) * symbol_length + 1;

    % gdzie się kończy
    end_index = start_index + symbol_length - 1;

    % wycinamy symbol
    segment = signal(start_index : end_index);

    % Obliczamy korelacje segmentu z każdym wzorcem
    corr0 = corr(segment.', wave0.');
    corr90 = corr(segment.', wave90.');
    corr180 = corr(segment.', wave180.');
    corr270 = corr(segment.', wave270.');

    % Szukamy największej korelacji i wybieramy indeks
    [~, idx] = max([corr0, corr90, corr180, corr270]);

    switch idx
        case 1
            recovered_bits = [recovered_bits '00'];
        case 2
            recovered_bits = [recovered_bits '01'];
        case 3
            recovered_bits = [recovered_bits '11'];
        case 4
            recovered_bits = [recovered_bits '10'];
    end
end

% Zamiana bitów na znaki
recovered_chars = char(bin2dec(reshape(recovered_bits,8,[])')).';
disp(['Odzyskane znaki: ' recovered_chars]);
