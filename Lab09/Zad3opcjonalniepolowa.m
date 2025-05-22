%% Cyfrowa PLL dla pilota 19 kHz
clear; close all; clc;

%% Parametry
fs = 48000;             % częstotliwość próbkowania, Hz
fpilot = 19000;         % częstotliwość pilota, Hz
omega = 2*pi*fpilot/fs; % początkowa estymacja częstotliwości kątowej (rad/próbkę)
alpha = 1e-2;           % współczynnik filtru IIR (faza)
beta = alpha^2/4;       % współczynnik filtru IIR (częstotliwość)
dur = 10;               % czas trwania sygnału, s
t = (0:1/fs:dur-1/fs)'; % wektor czasu
tolHz = 0.1;            % dokładność synchronizacji, Hz
consCnt = fs*0.01;      % liczba próbek dla stabilnej synchronizacji (0.01 s)

%% Funkcja uruchamiająca PLL
function [fEst, thetaEst, lockIdx, phase_error] = runPLL(p, fs, omega, alpha, beta, fpilot, tolHz, consCnt)
    N = length(p);
    theta = zeros(N+1, 1);
    freq = omega;
    fEst = zeros(N, 1);
    phase_error = zeros(N, 1);
    lockIdx = NaN;
    cntOK = 0;
    
    for n = 1:N
        perr = -p(n)*sin(theta(n));  % detektor fazy
        phase_error(n) = perr;        % zapisujemy błąd fazowy do analizy
        theta(n+1) = theta(n) + freq + alpha*perr;  % aktualizacja fazy
        freq = freq + beta*perr;      % aktualizacja częstotliwości
        fEst(n) = freq*fs/(2*pi);     % bieżąca estymacja w Hz
        
        % sprawdzanie synchronizacji
        if abs(fEst(n) - fpilot) < tolHz
            cntOK = cntOK + 1;
            if isnan(lockIdx) && cntOK >= consCnt
                lockIdx = n - consCnt + 1;
            end
        else
            cntOK = 0;
        end
    end
    
    thetaEst = theta(1:end-1);
end

%% 1) Syntetyczny pilot 19 kHz z przesunięciem fazowym
phi0 = 1.234;  % dowolne przesunięcie fazowe
p1 = cos(2*pi*fpilot*t + phi0);  % sygnał pilota
[f1, th1, lock1, perr1] = runPLL(p1, fs, omega, alpha, beta, fpilot, tolHz, consCnt);
fprintf('1) Synchronizacja na syntetycznym pilocie: %d próbek (%.3f s)\n', lock1, lock1/fs);

figure('Name', 'PLL - Pilot o stałej częstotliwości');

% Estymacja częstotliwości
subplot(2, 1, 1);
plot(t, f1, 'b'); hold on;
plot(t, fpilot*ones(size(t)), 'r--');
title('PLL na syntetycznym pilocie 19 kHz');
xlabel('Czas [s]');
ylabel('Estymowana częstotliwość [Hz]');
legend('f_{est}', 'f_{pilot}');
grid on;

% Porównanie faz
subplot(2, 1, 2);
ref_phase = mod(2*pi*fpilot*t + phi0, 2*pi);
pll_phase = mod(th1, 2*pi);
phase_diff = unwrap(pll_phase - ref_phase);
plot(t, phase_diff);
title('Różnica fazowa między sygnałem referencyjnym a odtworzonym');
xlabel('Czas [s]');
ylabel('Różnica fazy [rad]');
grid on;

%% 2) Pilot z wolną modulacją ±10 Hz, fm=0.1 Hz
df = 10;  % dewiacja, Hz
fm = 0.1;  % częstotliwość modulacji, Hz
fmod = fpilot + df*sin(2*pi*fm*t);  % częstotliwość modulowana
p2 = cos(2*pi*cumtrapz(t, fmod) + phi0);  % całkowanie dla fazy
[f2, th2, lock2, perr2] = runPLL(p2, fs, omega, alpha, beta, fpilot, tolHz, consCnt);
fprintf('2) Synchronizacja przy wolnej FM: %d próbek (%.3f s)\n', lock2, lock2/fs);

figure('Name', 'PLL - Pilot o zmiennej częstotliwości');

% Porównanie częstotliwości
subplot(3, 1, 1);
plot(t, fmod, 'k', t, f2, 'b');
title('PLL przy FM \pm10 Hz, fm=0.1 Hz');
xlabel('Czas [s]');
ylabel('Częstotliwość [Hz]');
legend('f_{rzeczywista}', 'f_{estymowana}');
grid on;

% Błąd częstotliwości
subplot(3, 1, 2);
plot(t, f2 - fmod);
title('Błąd estymacji częstotliwości');
xlabel('Czas [s]');
ylabel('Błąd [Hz]');
grid on;

% Błąd fazowy
subplot(3, 1, 3);
plot(t, perr2);
title('Błąd fazowy w czasie');
xlabel('Czas [s]');
ylabel('Błąd fazowy');
grid on;

%% 3) Szybkość synchronizacji przy różnych SNR
SNRlist = [Inf, 20, 10, 5, 0];  % dB, zgodnie z wymaganiami zadania
lockSNR = NaN(size(SNRlist));
convergence_ms = NaN(size(SNRlist));

% Przygotowanie wykresu
figure('Name', 'PLL - Szybkość zbieżności');

% Pętla dla różnych wartości SNR
for k = 1:numel(SNRlist)
    snrVal = SNRlist(k);
    
    % Generacja zaszumionego sygnału
    if isfinite(snrVal)
        p_noisy = awgn(p1, snrVal, 'measured');  % dodanie szumu
    else
        p_noisy = p1;  % czyste bez szumu
    end
    
    % Uruchomienie PLL
    [f_snr, th_snr, lock_snr, perr_snr] = runPLL(p_noisy, fs, omega, alpha, beta, fpilot, tolHz, consCnt);
    lockSNR(k) = lock_snr;
    convergence_ms(k) = lock_snr/fs*1000;  % czas w ms
    
    % Wykresy dla każdego SNR
    subplot(length(SNRlist), 1, k);
    t_short = t(1:min(length(t), fs/2));  % pokazujemy pierwsze 0.5s
    
    % Indeksy dla skróconego zakresu czasu
    idx_short = 1:min(length(t), fs/2);
    
    plot(t_short, p_noisy(idx_short), 'b'); 
    hold on;
    plot(t_short, cos(th_snr(idx_short)), 'r--');
    
    % Zaznaczenie momentu synchronizacji jeśli wystąpił w widocznym zakresie
    if ~isnan(lock_snr) && lock_snr <= fs/2
        plot([t(lock_snr), t(lock_snr)], [-1.5, 1.5], 'g-', 'LineWidth', 2);
    end
    
    hold off;
    grid on;
    
    if isinf(snrVal)
        title_str = ['SNR = \infty dB, Czas zbieżności = ', num2str(convergence_ms(k), '%.2f'), ' ms'];
    else
        title_str = ['SNR = ', num2str(snrVal), ' dB, Czas zbieżności = ', num2str(convergence_ms(k), '%.2f'), ' ms'];
    end
    
    title(title_str);
    xlabel('Czas [s]');
    ylabel('Amplituda');
    ylim([-1.5, 1.5]);
    legend('Sygnał zaszumiony', 'Sygnał odtworzony przez PLL');
end

% Tabela czasów zbieżności
figure('Name', 'PLL - Czasy zbieżności');

% Wykres słupkowy
bar(convergence_ms);

% Etykiety osi X
x_labels = cell(size(SNRlist));
for i = 1:length(SNRlist)
    if isinf(SNRlist(i))
        x_labels{i} = '∞';
    else
        x_labels{i} = num2str(SNRlist(i));
    end
end

set(gca, 'XTickLabel', x_labels);
grid on;
title('Czasy zbieżności PLL dla różnych wartości SNR');
xlabel('SNR [dB]');
ylabel('Czas zbieżności [ms]');

% Tabela z wynikami
Tlock = lockSNR/fs;  % czas, s
convergence_ms = Tlock * 1000;  % czas, ms
Ttbl = table(SNRlist', lockSNR', Tlock', convergence_ms', ...
    'VariableNames', {'SNR_dB', 'Liczba_probek_synch', 'Czas_synch_s', 'Czas_synch_ms'});

disp('3) Czas synchronizacji przy różnych SNR:');
disp(Ttbl);



%% 4) Separacja stereo FM bez PLL – idealny c38 = cos(2*pi*38kHz*t)

% Zakładamy, że 'y' to zdemodulowany sygnał FM (zawiera L+R, pilot 19kHz, L-R)
% Użyj pasmowego filtru hBP38 wokół 38kHz i LPF hLP (do ~15kHz)

t = (0:length(y)-1)'/fs;        % wektor czasu
fpilot = 19000;
c38 = cos(2*pi*2*fpilot*t);     % idealny sygnał 38 kHz

% Ekstrakcja L-R (z pasma 38 kHz)
y_LmR = filter(hBP38, 1, y);
LmR = 2 * y_LmR .* c38;         % demodulacja DSB-SC
LmR = filter(hLP, 1, LmR);      % pasmo audio

% L+R z basebandu
LpR = filter(hLP, 1, y);

% Rekonstrukcja kanałów stereo
L = (LpR + LmR)/2;
R = (LpR - LmR)/2;

figure('Name', 'Punkt 4 – Separacja kanałów bez PLL');
subplot(2,1,1); plot(t, L); title('Kanał L (bez PLL)'); grid on;
subplot(2,1,2); plot(t, R); title('Kanał R (bez PLL)'); grid on;


%% 5) Separacja kanałów przy odstrojonym lub przesuniętym c38

% Przypadek A – odstrojenie c38 (np. 38001 Hz)
fpilot_off = 19000.5;                 % 38001 Hz
c38_off = cos(2*pi*2*fpilot_off*t);
LmR_off = 2 * filter(hBP38,1,y) .* c38_off;
LmR_off = filter(hLP, 1, LmR_off);
L_off = (LpR + LmR_off)/2;
R_off = (LpR - LmR_off)/2;

% Przypadek B – przesunięcie fazy (np. 60 stopni)
phi = pi/3;
c38_phi = cos(2*pi*2*fpilot*t + phi);
LmR_phi = 2 * filter(hBP38,1,y) .* c38_phi;
LmR_phi = filter(hLP, 1, LmR_phi);
L_phi = (LpR + LmR_phi)/2;
R_phi = (LpR - LmR_phi)/2;

% Wykres
figure('Name', 'Punkt 5 – Separacja przy błędach c38');
subplot(2,2,1); plot(t, L_off); title('L przy 38001 Hz'); grid on;
subplot(2,2,2); plot(t, R_off); title('R przy 38001 Hz'); grid on;
subplot(2,2,3); plot(t, L_phi); title('L przy φ = 60°'); grid on;
subplot(2,2,4); plot(t, R_phi); title('R przy φ = 60°'); grid on;


