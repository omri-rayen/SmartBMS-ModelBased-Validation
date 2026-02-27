clear; clc; close all;

%% Paramètres
N_series      = 90;
R0_correction = 0.015;              % Correction ohmique appliquée : V_ocv = V_cell_mesurée - I * R0_correction (Ohm)
seuil_courant = 0.1;                % (A)
Ts            = 8;                  % Période d'échantillonnage des données (s)
dwell_samples = ceil(600 / Ts);     % Nombre d'échantillons correspondant à 600s (10 min) de "dwell" requis pour relaxation
T_min         = 20;
T_max         = 30;
min_pts_bin   = 50;                 % Nombre minimal de points requis dans un bin de SoC pour accepter la médiane
num_files     = 20;

save_path = fullfile('..', 'data', 'parameters');
if ~exist(save_path, 'dir'), mkdir(save_path); end

% Référence CATL NMC
soc_nmc_ref = [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100]';
ocv_nmc_ref = [3.00, 3.18, 3.38, 3.52, 3.60, 3.66, 3.72, 3.79, 3.87, 3.97, 4.09, 4.16, 4.20]';

%% Extraction des points de relaxation
all_soc = []; all_ocv = [];   % valeurs SoC & OCV valides
files_found = 0;              % fichiers lus

for f = 1:num_files
    filename = fullfile('..', 'data', 'raw', sprintf('#%d.csv', f));
    if ~exist(filename, 'file'), continue; end

    data = readtable(filename);
    files_found = files_found + 1;

    data.Properties.VariableNames(1:11) = {'number','record_time','soc','pack_voltage',...
        'charge_current','max_cell_voltage','min_cell_voltage',...
        'max_temperature','min_temperature','available_energy','available_capacity'};

    cell_v  = data.pack_voltage / N_series;                             % tension par cellule (vecteur)
    low_I   = abs(data.charge_current) < seuil_courant;                 % 0 (ex: -5.0) ou 1 (ex: 0.03)
    % détecter une séquence continue de courant quasi nul pendant dwell_samples points
    relaxed = movmin(double(low_I), [dwell_samples - 1, 0]) == 1;       % movmin(signal, [Nb_avant, Nb_après])
    temp_ok = data.max_temperature >= T_min & data.max_temperature <= T_max;

    idx   = relaxed & temp_ok;                                          % index des points au repos
    ocv_c = cell_v(idx) - data.charge_current(idx) * R0_correction;

    all_soc = [all_soc; data.soc(idx)];
    all_ocv = [all_ocv; ocv_c];

    fprintf('#%d : %d points\n', f, sum(idx));
end
fprintf('%d fichiers, %d points au total\n', files_found, length(all_soc));

%% Médiane par bin + filtrage par rapport à la référence
soc_bins   = 0:1:100;                           % bins largeur 1%
soc_center = (soc_bins(1:end-1) + 0.5)';        % (ex: 0.5, 1.5, 2.5)
ocv_median = NaN(size(soc_center));
ocv_count  = zeros(size(soc_center));

for i = 1:length(soc_center)
    % indices des points expérimentaux tombant dans le bin i
    idx = all_soc >= soc_bins(i) & all_soc < soc_bins(i+1);
    if sum(idx) >= min_pts_bin
        m_val = median(all_ocv(idx));
        % récupèrer la valeur de référence (extrapolation linéaire si nécessaire)
        v_ref = interp1(soc_nmc_ref, ocv_nmc_ref, soc_center(i), 'linear', 'extrap');
        if abs(m_val - v_ref) < 0.12   % tolérance 120 mV
            ocv_median(i) = m_val;
            ocv_count(i) = sum(idx);
        end
    end
end

% Sélection des bins valides (pas NaN)
valid = ~isnan(ocv_median);
s_v = soc_center(valid);    % SoC validés (centres)
v_v = ocv_median(valid);    % OCV médians correspondants

%% Fusion données expérimentales + référence, puis lissage
soc_merge = [s_v; soc_nmc_ref];
ocv_merge = [v_v; ocv_nmc_ref];
[soc_merge, ia] = sort(soc_merge);
ocv_merge = ocv_merge(ia);
% Supprimer les duplicatas de SoC proches (arrondis à 0.1)
[soc_merge, ib] = unique(round(soc_merge, 1));
ocv_merge = ocv_merge(ib);

soc_interp = (0:1:100)';
% Interpole (pchip pour garder forme lisse) puis lissage local (loess sur fenêtre de 15 points)
ocv_interp = smoothdata(interp1(soc_merge, ocv_merge, soc_interp, 'pchip'), 'loess', 15);

% Monotonie stricte (requis pour EKF)
ocv_interp(1) = ocv_nmc_ref(1);
ocv_interp(end) = ocv_nmc_ref(end);
pente_min = 0.001; 
for k = 2:length(ocv_interp)
    if ocv_interp(k) < (ocv_interp(k-1) + pente_min)
        ocv_interp(k) = ocv_interp(k-1) + pente_min;
    end
end

%% Sauvegarde
nmc_params.soc         = soc_interp;
nmc_params.ocv         = ocv_interp;
nmc_params.date_export = datetime('now');

save(fullfile(save_path, 'nmc_params.mat'), 'nmc_params');
fprintf('Sauvegardé dans nmc_params.mat\n');

%% Affichage
figure('Position', [100 100 1200 500]);

subplot(1,2,1); hold on; grid on;
plot(soc_nmc_ref, ocv_nmc_ref, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Référence CATL');
scatter(s_v, v_v, 30, 'b', 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', 'Expérimental');
plot(soc_interp, ocv_interp, 'r-', 'LineWidth', 2, 'DisplayName', 'Modèle lissé');
xlabel('SoC (%)'); ylabel('OCV (V)');
title('Courbe OCV-SoC');
legend('Location', 'southeast');
xlim([0 100]); ylim([2.7 4.3]);

subplot(1,2,2);
bar(soc_center, ocv_count, 'FaceColor', [0 0.447 0.741], 'EdgeColor', 'none');
hold on;
yline(min_pts_bin, 'r--', 'Seuil de validité', 'LineWidth', 1.5);
xlabel('SoC (%)'); ylabel('Nombre de points');
title('Points valides par bin');
grid on;