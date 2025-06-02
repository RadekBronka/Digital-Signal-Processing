function dq = lab11_kwant(d)
    % Kwantyzacja do 16 poziomów (4 bity)
    min_d = min(d);
    max_d = max(d);
    L = 256; % liczba poziomów
    step = (max_d - min_d) / (L - 1);  % szerokość przedziału

    % Kwantyzacja – zaokrąglenie do najbliższego poziomu
    dq = round((d - min_d) / step) * step + min_d;
end
