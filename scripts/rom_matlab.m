% rom_matlab.m - thermisches ROM des Schaltschranks in MATLAB (ohne Toolboxen).
% Liest ergebnisse/rom/rom_parameter.json (erzeugt von rom.py) und rechnet
%   (1) das nichtlineare RC-Netzwerk mit ode15s (Sprung aller Verluste aus dem kalten Zustand),
%   (2) das um den Basis-Betriebspunkt linearisierte Zustandsraummodell (expm, exakt diskretisiert).
% Ergebnis: ergebnisse/rom/rom_matlab.csv (t, 6 Knoten nichtlinear, 6 Knoten linear), Laufzeiten.
% Aufruf:  matlab -batch "rom_matlab"   (im Projektordner)

par = jsondecode(fileread(fullfile('ergebnisse', 'rom', 'rom_parameter.json')));
sigma = par.sigma;
Cap = par.C(:);
P = par.P_basis(:);
TU = par.T_U_basis;
seiten = fieldnames(par.A_seite);
hA = 0;
for i = 1:numel(seiten)
    hA = hA + par.h_basis.(seiten{i}) * par.A_seite.(seiten{i});
end
gL = par.gL(:); gW = par.gW(:); gLW = par.gLW;
f = @(t, T) rhs(T, P, TU, gL, gW, gLW, hA, par.eps_aussen, sigma, par.A_ab) ./ Cap;

t_end = 8 * 3600;
tt = (0:60:t_end)';
tic;
opt = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
[~, Tnl] = ode15s(f, tt, TU * ones(6, 1), opt);
t_nl = toc;

% linearisiert um den stationaeren Basispunkt: x = T - Tss, Eingangssprung ab kaltem Zustand
A = par.A_lin; B = par.B_lin; Tss = par.T_ss_basis(:);
tic;
dt = 60; Ad = expm(A * dt);
x = (TU - Tss);            % Start: alle Knoten auf Umgebungstemperatur (Abweichung vom Betriebspunkt)
x = x .* ones(6, 1);
Tlin = zeros(numel(tt), 6);
for k = 1:numel(tt)
    Tlin(k, :) = (Tss + x)';
    x = Ad * x;            % Eingaenge = Betriebspunkt -> nur Eigenbewegung um Tss
end
t_lin = toc;

M = [tt, Tnl - 273.15, Tlin - 273.15];
kopf = 't_s,U_nl,N_nl,S_nl,T_nl,Luft_nl,Wand_nl,U_lin,N_lin,S_lin,T_lin,Luft_lin,Wand_lin';
fid = fopen(fullfile('ergebnisse', 'rom', 'rom_matlab.csv'), 'w');
fprintf(fid, '%s\n', kopf);
fclose(fid);
writematrix(M, fullfile('ergebnisse', 'rom', 'rom_matlab.csv'), 'WriteMode', 'append');
fprintf('MATLAB-ROM: nichtlinear (ode15s) %.3f s, linear (expm) %.3f s fuer %.0f h\n', t_nl, t_lin, t_end / 3600);
fprintf('Endwerte nichtlinear [degC]: %s\n', sprintf('%.2f ', Tnl(end, :) - 273.15));
fprintf('Endwerte linear      [degC]: %s\n', sprintf('%.2f ', Tlin(end, :) - 273.15));

function r = rhs(T, P, TU, gL, gW, gLW, hA, epsA, sigma, Aab)
    Ti = T(1:4); TL = T(5); TW = T(6);
    dL = Ti - TL;
    qiL = gL .* (dL .^ 2 + 1) .^ (1/6) .* dL;
    qiW = gW .* sigma .* (Ti .^ 4 - TW ^ 4);
    dLW = TL - TW;
    qLW = gLW * (dLW ^ 2 + 1) ^ (1/6) * dLW;
    qWU = hA * (TW - TU) + epsA * sigma * Aab * (TW ^ 4 - TU ^ 4);
    r = zeros(6, 1);
    r(1:4) = P - qiL - qiW;
    r(5) = sum(qiL) - qLW;
    r(6) = sum(qiW) + qLW - qWU;
end
