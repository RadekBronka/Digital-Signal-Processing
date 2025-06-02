[x, fs] = audioread('DontWorryBeHappy.wav', 'native'); % wczytanie próbki
x = double(x); % zamiana na double
% x jest teraz Nx2 (stereo)

a = 0.9545; % parametr kodera

% przygotowanie wyjścia
y = zeros(size(x));
dq = zeros(size(x));

for ch = 1:2
    % --- KODER ---
    d = x(:,ch) - a * [0; x(1:end-1,ch)];
    
    % --- KWANTYZATOR ---
    dq(:,ch) = lab11_kwant(d);
    
    % --- DEKODER ---
    y(:,ch) = zeros(size(x,1),1);
    for n = 2:length(x)
        y(n,ch) = dq(n,ch) + a * y(n-1,ch);
    end
end

% Przycięcie do 5 sekund
samples5s = min(length(x), fs * 5);

x = x(1:samples5s, :);
dq = dq(1:samples5s, :);
y = y(1:samples5s, :);

n = 1:samples5s;

% --- WYKRES 1: Oryginał vs dq (przed dekoderem) dla kanału lewego ---
figure(1);
plot(n, x(:,1), 'b', n, dq(:,1), 'r');
legend('x(n) lewy - oryginał', 'dq(n) lewy - przed dekoderem');
title('Oryginał vs sygnał przed dekoderem (kanał lewy)');
xlabel('Numer próbki');
ylabel('Amplituda');

% --- WYKRES 2: Oryginał vs y (po dekodowaniu) dla kanału lewego ---
figure(2);
plot(n, x(:,1), 'b', n, y(:,1), 'g');
legend('x(n) lewy - oryginał', 'y(n) lewy - po dekodowaniu');
title('Oryginał vs sygnał po dekodowaniu (kanał lewy)');
xlabel('Numer próbki');
ylabel('Amplituda');

% --- WYKRES 3: Oryginał vs y (po dekodowaniu) dla kanału prawego ---
figure(3);
plot(n, x(:,2), 'b', n, y(:,2), 'g');
legend('x(n) prawy - oryginał', 'y(n) prawy - po dekodowaniu');
title('Oryginał vs sygnał po dekodowaniu (kanał prawy)');
xlabel('Numer próbki');
ylabel('Amplituda');

% Odtwarzanie dźwięku stereo po przetworzeniu
sound(y, fs);
