fpr=16000;
T=0.1;
fc=500;

dt=1/fpr;

t=0:dt:T-dt;

signal=[];

imie='Radek';
imied=double(imie);
imiebin=dec2bin(imied);

ciag = reshape(imiebin.', 1, []); %zamieniam macierz imiebin na ciag 
% (1- jeden wiersz, []- tyle kolunm ile potrzeba- oblicza matlab

for bit = ciag
    if bit == '0'
        signal = [signal, sin(2 * pi * fc * t)];
    else
        signal = [signal, -sin(2 * pi * fc * t)];
    end
end

% Rysowanie sygnału
plot(signal);
title('Sygnał transmisji bitów imienia Radek');
xlabel('Próbka');
ylabel('Amplituda');

% Odsłuchanie sygnału w różnych częstotliwościach próbkowania
soundsc(signal, 8000);
pause(1.5);
soundsc(signal, 16000);
pause(1.5);
soundsc(signal, 24000);
pause(1.5);
soundsc(signal, 32000);
pause(1.5);
soundsc(signal, 48000);


