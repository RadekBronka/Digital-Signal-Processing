clear all;
close all;
clearvars;
clc;

load("adsl_x.mat");
M=32;
N=512;
K=4;

positions=[];

for i = N:N:length(x)
    last32 = x(i-N+1:i); %ostatnie 32 elementy bloku
    %korelacja
    [corr_vals, lags] = xcorr(x, last32);
    
    [~, idx] = max((corr_vals));  % Indeks maksymalnej korelacji

    % Przechowywanie pozycji początkowej prefiksu
    position = lags(idx)+1;
    positions = [positions, position];

end
