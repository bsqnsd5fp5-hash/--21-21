% Этап 2: Моделирование траектории спутника GPS PRN 24
% Эфемериды взяты из out.txt (Этап 1)
clear; clc; close all;

% ----- Эфемериды спутника 24 -----
ephem.Crs      = 2.275000e+01;
ephem.Dn       = 5.210052e-07;     % град/с
ephem.M0       = -12.755586;       % град
ephem.Cuc      = -1.557082e-06;
ephem.e        = 1.636439e-02;
ephem.Cus      = 6.183982e-06;
ephem.sqrtA    = 5.153701e+03;
ephem.toe      = 288000;           % секунды недели GPS
ephem.Cic      = -2.142042e-07;
ephem.Omega0   = -147.307115;      % град
ephem.Cis      = 1.266599e-07;
ephem.i0       = 53.402052;        % град
ephem.Crc      = 1.897500e+02;
ephem.omega    = 82.744271;        % град
ephem.OmegaDot = -4.756403e-07;    % град/с
ephem.iDot     = -4.752011e-08;    % град/с
ephem.toc      = 288000;
ephem.af2      = 0.0;
ephem.af1      = 2.159915e-12;
ephem.af0      = -2.265218e-04;

% ----- Константы -----
GM   = 3.986005e14;              % гравитационная постоянная Земли
OMGE = 7.2921151467e-5;          % угловая скорость Земли, рад/с
DEG2RAD = pi/180;
R_EARTH = 6371000;               % средний радиус Земли, м

% ----- Параметры наблюдателя (Москва, МЭИ) -----
lat = 55.7540 * DEG2RAD;
lon = 37.6200 * DEG2RAD;
h_obs = 200;                     % высота над эллипсоидом, м

% ----- Интервал моделирования -----
% 3 марта 2021 г. (день недели 3), начало дня в секундах недели = 3*86400
start_sec = 259200;
end_sec   = 345600;
dt = 300;                        % шаг 5 минут
t_week = start_sec:dt:end_sec;   % секунды недели

% ----- Расчёт всех положений -----
N = length(t_week);
x_ecef = zeros(1, N);
y_ecef = zeros(1, N);
z_ecef = zeros(1, N);

for i = 1:N
    [x_ecef(i), y_ecef(i), z_ecef(i)] = satpos(ephem, t_week(i));
end

% ----- Запись координат в файл out_xyz.txt -----
sec_of_day = (t_week - start_sec);   % секунды от начала суток
fid = fopen('out_xyz.txt', 'w');
for i = 1:N
    fprintf(fid, '%d %.3f %.3f %.3f\n', sec_of_day(i), x_ecef(i), y_ecef(i), z_ecef(i));
end
fclose(fid);
disp('out_xyz.txt записан.');

% ----- 3D-график в ECEF -----
figure('Color','white');
plot3(x_ecef, y_ecef, z_ecef, 'b-', 'LineWidth', 1);
hold on;
% Сфера Земли
[u, v] = meshgrid(linspace(0, 2*pi, 60), linspace(0, pi, 60));
xs = R_EARTH * cos(u) .* sin(v);
ys = R_EARTH * sin(u) .* sin(v);
zs = R_EARTH * ones(size(u)) .* cos(v);
surf(xs, ys, zs, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', 'g');
xlabel('X (м)');
ylabel('Y (м)');
zlabel('Z (м)');
title('Траектория спутника PRN 24 в ECEF (3.03.2021)');
legend({'PRN 24', 'Земля'}, 'Location', 'best');
axis equal;
grid on;
saveas(gcf, 'ecef_trajectory.png');
disp('График ECEF сохранён (ecef_trajectory.png).');

% ----- ECI (инерциальная СК) -----
theta = OMGE * (t_week - start_sec);   % угол поворота Земли от начала дня
cos_th = cos(theta);
sin_th = sin(theta);
x_eci = x_ecef .* cos_th - y_ecef .* sin_th;
y_eci = x_ecef .* sin_th + y_ecef .* cos_th;
z_eci = z_ecef;

figure('Color','white');
plot3(x_eci, y_eci, z_eci, 'r-', 'LineWidth', 1);
hold on;
surf(xs, ys, zs, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', 'g');
xlabel('X (м)');
ylabel('Y (м)');
zlabel('Z (м)');
title('Траектория спутника PRN 24 в ECI (3.03.2021)');
legend({'PRN 24', 'Земля'}, 'Location', 'best');
axis equal;
grid on;
saveas(gcf, 'eci_trajectory.png');
disp('График ECI сохранён (eci_trajectory.png).');

% ----- SkyPlot -----
az = zeros(1, N);
el = zeros(1, N);
rx0 = R_EARTH * cos(lat) * cos(lon);
ry0 = R_EARTH * cos(lat) * sin(lon);
rz0 = R_EARTH * sin(lat);

for i = 1:N
    dx = x_ecef(i) - rx0;
    dy = y_ecef(i) - ry0;
    dz = z_ecef(i) - rz0;
    e = -sin(lon)*dx + cos(lon)*dy;
    n = -sin(lat)*cos(lon)*dx - sin(lat)*sin(lon)*dy + cos(lat)*dz;
    u = cos(lat)*cos(lon)*dx + cos(lat)*sin(lon)*dy + sin(lat)*dz;
    az(i) = atan2d(e, n);    % азимут в градусах
    az(i) = mod(az(i), 360);
    el(i) = asind(u / sqrt(e^2 + n^2 + u^2));  % угол места
end

figure('Color','white');
polarplot(deg2rad(az), 90 - el, 'b-', 'LineWidth', 1.5);
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.RLim = [0 90];
ax.RTick = [30 60 90];
% Убираем числовые метки у радиальных делений (если нужно)
set(ax, 'RTickLabel', []);
title('SkyPlot PRN 24 (3.03.2021, Москва)');
legend('PRN 24', 'Location', 'best');
grid on;
saveas(gcf, 'skyplot_matlab.png');
disp('SkyPlot сохранён (skyplot_matlab.png).');

disp('Все графики построены.');

% ----- Функция расчёта положения спутника -----
function [xECEF, yECEF, zECEF] = satpos(ephem, t)
    GM = 3.986005e14;
    OMGE = 7.2921151467e-5;
    DEG2RAD = pi/180;
    
    tk = t - ephem.toe;
    if tk > 302400
        tk = tk - 604800;
    elseif tk < -302400
        tk = tk + 604800;
    end
    n0 = sqrt(GM / ephem.sqrtA^6);
    n = n0 + ephem.Dn * DEG2RAD;
    M = ephem.M0 * DEG2RAD + n * tk;
    % Уравнение Кеплера
    E = M;
    for j = 1:10
        E_old = E;
        E = M + ephem.e * sin(E_old);
        if abs(E - E_old) < 1e-12
            break;
        end
    end
    v = 2 * atan2(sqrt(1+ephem.e)*sin(E/2), sqrt(1-ephem.e)*cos(E/2));
    phi = v + ephem.omega * DEG2RAD;
    du = ephem.Cus*sin(2*phi) + ephem.Cuc*cos(2*phi);
    dr = ephem.Crs*sin(2*phi) + ephem.Crc*cos(2*phi);
    di = ephem.Cis*sin(2*phi) + ephem.Cic*cos(2*phi);
    u = phi + du;
    r = ephem.sqrtA^2 * (1 - ephem.e*cos(E)) + dr;
    i_orbit = ephem.i0*DEG2RAD + di + ephem.iDot*DEG2RAD*tk;
    x_orb = r * cos(u);
    y_orb = r * sin(u);
    OMEGA = ephem.Omega0*DEG2RAD + (ephem.OmegaDot*DEG2RAD - OMGE)*tk - OMGE*ephem.toe;
    xECEF = x_orb*cos(OMEGA) - y_orb*cos(i_orbit)*sin(OMEGA);
    yECEF = x_orb*sin(OMEGA) + y_orb*cos(i_orbit)*cos(OMEGA);
    zECEF = y_orb*sin(i_orbit);
end