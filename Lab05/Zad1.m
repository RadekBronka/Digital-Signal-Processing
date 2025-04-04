%% Parametry zer i biegunów
% Zera: z = ±j5, ±j15
zeros_tf = [1j*5, -1j*5, 1j*15, -1j*15];

% Bieguny: p1,2 = -0.5 ± j9.5; p3,4 = -1 ± j10; p5,6 = -0.5 ± j10.5
poles_tf = [ -0.5+1j*9.5, -0.5-1j*9.5, -1+1j*10, -1-1j*10, -0.5+1j*10.5, -0.5-1j*10.5];

%% Wyznaczenie wielomianów licznikowego i mianownikowego
num_poly = poly(zeros_tf);   % współczynniki licznika
den_poly = poly(poles_tf);     % współczynniki mianownika

%% Dobór wzmocnienia k
% Chcemy uzyskać |H(j10)| = 1, czyli przy s = j*10
w_pass = 10;
s = 1j*w_pass;
H_j10 = polyval(num_poly, s) / polyval(den_poly, s);
k = 1/abs(H_j10);

% Wyznaczenie ostatecznych współczynników transmitancji
num_coeff = k * num_poly;
den_coeff = den_poly;

fprintf('Dobór k = %f\n', k);

%% Definicja funkcji transmitancji H(jw)
% Obliczamy H(jw) dla wektora częstotliwości
omega = linspace(0,40,1000);
s = 1j*omega;
H = polyval(num_coeff, s) ./ polyval(den_coeff, s);

%% Charakterystyka amplitudowo-częstotliwościowa
amplitude = abs(H);
amplitude_dB = 20*log10(amplitude);

figure;
subplot(2,1,1);
plot(omega, amplitude, 'LineWidth',1.5);
grid on;
xlabel('\omega [rad/s]');
ylabel('|H(j\omega)|');
title('Charakterystyka amplitudowa (skala liniowa)');

subplot(2,1,2);
plot(omega, amplitude_dB, 'LineWidth',1.5);
grid on;
xlabel('\omega [rad/s]');
ylabel('20 log|H(j\omega)| [dB]');
title('Charakterystyka amplitudowa (skala decybelowa)');

%% Charakterystyka fazowo-częstotliwościowa
phase_deg = unwrap(angle(H))*180/pi; % faza w stopniach
figure;
plot(omega, phase_deg, 'LineWidth',1.5);
grid on;
xlabel('\omega [rad/s]');
ylabel('Faza [°]');
title('Charakterystyka fazowa');

%% Rozmieszczenie biegunów i zer na płaszczyźnie zespolonej
figure;
plot(real(zeros_tf), imag(zeros_tf), 'o', 'MarkerSize',10, 'LineWidth',2);
hold on;
plot(real(poles_tf), imag(poles_tf), '*', 'MarkerSize',10, 'LineWidth',2);
xlabel('Re');
ylabel('Im');
title('Rozmieszczenie biegunów i zer');
legend('Zera','Bieguny');
grid on;
axis equal;

%% Wyznaczenie tłumienia w paśmie zaporowym
% Zakładamy, że pasmo zaporowe to: omega < 5 rad/s oraz omega > 15 rad/s
idx_low = find(omega < 5);
idx_high = find(omega > 15);

% Tłumienie w dB
amp_dB_low = amplitude_dB(idx_low);
amp_dB_high = amplitude_dB(idx_high);

% W każdej z tych grup:
% Minimalne tłumienie - wartość najbliższa 0 dB (najmniejsze tłumienie)
% Maksymalne tłumienie - wartość najbardziej ujemna (największe tłumienie)
min_att_low = max(amp_dB_low); % najmniejsze tłumienie (mniej tłumione) w paśmie niskich częstotliwości
max_att_low = min(amp_dB_low); % największe tłumienie (bardziej tłumione) w paśmie niskich częstotliwości

min_att_high = max(amp_dB_high);
max_att_high = min(amp_dB_high);

fprintf('\nW paśmie zaporowym:\n');
fprintf('Dla omega < 5 rad/s: Minimalne tłumienie = %f dB, Maksymalne tłumienie = %f dB\n', min_att_low, max_att_low);
fprintf('Dla omega > 15 rad/s: Minimalne tłumienie = %f dB, Maksymalne tłumienie = %f dB\n', min_att_high, max_att_high);

% Ogólne wartości (można je potraktować jako najłagodniejsze i najbardziej tłumione miejsce w paśmie zaporowym)
overall_min_att = max([min_att_low, min_att_high]); % najmniejsze tłumienie w paśmie zaporowym (najmniej tłumione)
overall_max_att = min([max_att_low, max_att_high]); % największe tłumienie w paśmie zaporowym (najbardziej tłumione)

fprintf('Ogólnie w paśmie zaporowym: Minimalne tłumienie = %f dB, Maksymalne tłumienie = %f dB\n', overall_min_att, overall_max_att);
