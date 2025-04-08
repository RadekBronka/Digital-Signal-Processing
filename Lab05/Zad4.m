% cps_07_analog_tranform.m
clear all; close all;
% Wymagania
N = 8; % liczba biegunow transmitanci prototypu analogowego
f0 = 100; % dla filtrow LowPass oraz HighPass
f1 = 10; f2=100; % dla filtrow BandPass oraz BandStop
Rp = 3; % dozwolony poziom oscylacji w pasmie przepustowym (w dB) - ripples in pass
Rs = 100; % dozwolony poziom oscylacji w pasmie zaporowym (w dB) - ripples in stop

% Projekt analogowego filtra prototypowego - dolnoprzepustowy z w0 = 1
%[z,p,wzm] = buttap(N); % analogowy prototyp Butterwotha
 [z,p,wzm] = cheb1ap(N,Rp); % analogowy prototyp Czebyszewa typu-I
% [z,p,wzm] = cheb2ap(N,Rs); % analogowy prototyp Czebyszewa typu-II
% [z,p,wzm] = ellipap(N,Rp,Rs); % analogowy prototyp Cauera (eliptyczny)
b = wzm*poly(z); % [z,wzm] --> b,  poly zwraca współczynniki wielomianu
a = poly(p); % p --> a
f = 0 : 0.01 : 1000; w=2*pi*f; % zakres czestotliwosci
H = freqs(b,a,w); % odpowiedz czestotliwosciowa filtra (funkcja Matlaba) H(jw)
fi=0:pi/1000:2*pi; c=cos(fi); s=sin(fi);
figure; semilogx(w,20*log10(abs(H))); grid; xlabel('w [rad*Hz]'); title('Analog Proto |H(f)');
figure; plot(real(z),imag(z),'ro',real(p),imag(p),'b*',c,s,'k-'); grid; title('Analog Proto ZP'); 

% Transformacja czestotliwosci: znormalizowany (w0=1) LowPass --> xxxPass, xxxStop
% Funkcje xx2yy z Signal Processing Toolbox
[b,a] = lp2lp(b,a,2*pi*f0); % LowPass na LowPass
% [b,a] = lp2hp(b,a,2*pi*f0); % LowPass na HighPass
% [b,a] = lp2bp(b,a,2*pi*sqrt(f1*f2),2*pi*(f2-f1)); % LowPass na BandPass
% [b,a] = lp2bs(b,a,2*pi*sqrt(f1*f2),2*pi*(f2-f1)); % LowPass na BandStop
z=roots(b); p=roots(a); %nowe zera i bieguny
% ... kontynuuj cps_07_analog_intro.m
% Zaprojektuj/dobierz wspolczynniki transmitancji filtra analogowego
if(0) % dobor wartosci wspolczynnikow wielomianow zmiennej 's' transmitancji
    b = [3,2]; % [ b1, b0 ]
    a = [4,3,2,1]; % [ a3, a2, a1, a0=1]
    z = roots(b); p = roots(a); % [b,a] --> [z,p]
else % dobor wartosci pierwiastkow wielomianow zmiennej 's' transmitancji, tutaj najpierw zera i bieguny
    wzm = 0.001;
    z = j*2*pi*[ 600,800 ]; z = [z conj(z)];
    p = [-1,-1] + j*2*pi*[100,200]; p = [p conj(p)];
    b = wzm*poly(z); a = poly(p); % [z,p] --> [b,a]
end
figure; plot(real(z),imag(z),'bo', real(p),imag(p),'r*'); grid;

title('Zera (o) i Bieguny (*)'); xlabel('Real()'); ylabel('Imag()'); 

% Weryfikacja odpowiedzi (charakterystyki) filtra: amplitudowej fazowej, impulsowej, skokowej
f = 0 : 0.1 : 1000; % czestotliwosc w hercach
w = 2*pi*f; % czestotliwosc katowa, pulsacja
s = j*w; % zmienna transformacji Laplace'a
H = polyval(b,s) ./ polyval(a,s); % h(s) --> H(f) dla s=j*w: iloraz dwoch wielomianow
figure; plot(f,20*log10(abs(H))); xlabel('f [Hz]'); title('|H(f)| [dB]'); grid; %amplitudowa
figure; plot(f,unwrap(angle(H))); xlabel('f [Hz]'); title('angle(H(f)) [rad]'); grid; %fazowa
figure; impulse(b,a); % odpowiedz filtra na pobudzenie impulsowe (z Control Toolbox)
figure; step(b,a); % odpowiedz filtra na pobudzenie skokowe (z Control Toolbox)
