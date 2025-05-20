clear;
clc;

% Parametry sygnału
fs = 228000;           % częstotliwość próbkowania [Hz]
fpilot = 19000;        % częstotliwość pilota [Hz]
f_offset = 0;          % ewentualne rozstrojenie [Hz]
phi_offset = pi/4;     % przesunięcie fazowe [rad]
t = 0:1/fs:0.1;        % wektor czasu (0.1 s)
p = cos(2*pi*(fpilot + f_offset)*t + phi_offset);  % Sygnał pilota

% Inicjalizacja PLL
freq = 2*pi*fpilot/fs;          % częstotliwość początkowa w rad/sample
theta = zeros(1, length(p)+1);  % wektor fazy
alpha = 1e-2;                   % współczynnik proporcjonalny
beta = alpha^2/4;               % współczynnik całkujący

% Pętla PLL
for n = 1:length(p)
    perr = -p(n)*sin(theta(n));
    theta(n+1) = theta(n) + freq + alpha*perr;
    freq = freq + beta*perr;
end

% Generacja trzeciej harmonicznej (np. 57 kHz)
c57 = cos(3 * theta(1:end-1));

% Oblicz błąd fazy
true_phase = 2*pi*(fpilot + f_offset)*t + phi_offset;
estimated_phase = theta(1:end-1);
phase_error = wrapToPi(true_phase - estimated_phase);

% Wykres błędu fazy
figure;
plot(t, phase_error);
xlabel('Czas [s]');
ylabel('Błąd fazy [rad]');
title('Błąd fazy między sygnałem a estymowaną fazą PLL');
grid on;

% Estymowana częstotliwość
inst_freq = diff(theta) * fs / (2*pi);

% Wykres estymowanej częstotliwości
figure;
plot(t, inst_freq);
yline(fpilot, 'r--', '19 kHz');
xlabel('Czas [s]');
ylabel('Częstotliwość [Hz]');
title('Estymowana częstotliwość przez PLL');
grid on;

% Parametry wykrywania stabilizacji
threshold = 0.01;                  % maksymalny błąd fazy [rad]
stable_time = 0.01;                % wymagany czas stabilności [s]
stable_samples = round(stable_time * fs);  % ilość próbek

% Szukanie momentu dostrojenia
for n = 1:(length(phase_error) - stable_samples)
    window = abs(phase_error(n:n + stable_samples - 1));
    if all(window < threshold)
        lock_sample = n;
        lock_time = t(lock_sample);
        fprintf('PLL dostroiło się po %.4f s (%d próbek).\n', lock_time, lock_sample);
        break;
    end
end

% Zaznaczenie momentu dostrojenia na wykresie błędu fazy
hold on;
xline(lock_time, 'r--', 'Moment dostrojenia');



%% Pilot z modulowaną częstotliwością (±10 Hz, 0.1 Hz)
df = 10;                 % maksymalne odchylenie częstotliwości [Hz]
fm = 0.1;                % częstotliwość modulacji [Hz]
f_inst = fpilot + df * sin(2*pi*fm*t);  % chwilowa częstotliwość
phi_mod = 2*pi*cumsum(f_inst)/fs + phi_offset;  % zintegrowana faza
p_mod = cos(phi_mod);    % sygnał pilota z modulowaną częstotliwością

% Inicjalizacja PLL
freq = 2*pi*fpilot/fs;           % początkowa estymacja częstotliwości
theta = zeros(1, length(p_mod)+1);
alpha = 1e-2;
beta = alpha^2/4;

% Pętla PLL dla sygnału modulowanego
for n = 1:length(p_mod)
    perr = -p_mod(n)*sin(theta(n));
    theta(n+1) = theta(n) + freq + alpha*perr;
    freq = freq + beta*perr;
end

% Estymowana faza i częstotliwość
estimated_phase = theta(1:end-1);
true_phase_mod = phi_mod;
phase_error_mod = wrapToPi(true_phase_mod - estimated_phase);
inst_freq_mod = diff(theta) * fs / (2*pi);

% Wykres błędu fazy
figure;
plot(t, phase_error_mod);
xlabel('Czas [s]');
ylabel('Błąd fazy [rad]');
title('Błąd fazy PLL przy sinusoidalnie modulowanej częstotliwości pilota');
grid on;

% Wykres estymowanej częstotliwości
figure;
plot(t, f_inst, 'k--', 'DisplayName', 'Prawdziwa f pilota');
hold on;
plot(t(1:length(inst_freq_mod)), inst_freq_mod, 'b', 'DisplayName', 'PLL - estymowana');
xlabel('Czas [s]');
ylabel('Częstotliwość [Hz]');
title('PLL śledzi wolno zmieniającą się częstotliwość pilota');
legend;
grid on;

% Parametry wykrywania stabilizacji dla modulowanego sygnału
threshold_mod = 0.01;                   % próg błędu fazy [rad]
stable_time_mod = 0.01;                 % wymagany czas stabilności [s]
stable_samples_mod = round(stable_time_mod * fs);  % ilość próbek

% Szukanie momentu dostrojenia w przypadku modulacji
for n = 1:(length(phase_error_mod) - stable_samples_mod)
    window_mod = abs(phase_error_mod(n:n + stable_samples_mod - 1));
    if all(window_mod < threshold_mod)
        lock_sample_mod = n;
        lock_time_mod = t(lock_sample_mod);
        fprintf('PLL (modulowany sygnał) dostroiło się po %.4f s (%d iteracji).\n', ...
                lock_time_mod, lock_sample_mod);
        break;
    end
end

% Zaznaczenie momentu dostrojenia na wykresie błędu fazy
figure(gcf);  % aktywuj ostatnią figurę
hold on;
xline(lock_time_mod, 'r--', 'Moment dostrojenia');
