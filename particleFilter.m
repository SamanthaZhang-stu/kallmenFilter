%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SC23219041 张彦 2024/4/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 清除工作空间、命令窗口和图形窗口
clear, clc, close all

% 定义船只位置的取值范围
xmin = -10;
xmax = 10;
% 定义船只速度的取值范围
xpmin = -0.3;
xpmax = 0.3;
% 粒子数
N = 10000;
% 系统参数
m = 1;          % 质量（kg）
kk = 1;         % 弹簧常数（N/m）
c = 0.5;        % 阻尼系数（N/s）
F0 = 10;        % 扰动力的最大幅值（N）
dT = 0.05;      % 时间步长（s）
M = 400;        % 时间步数

% 计算真实的扰动力序列
wk = (rand(M, 1) - 0.5) * 2 * F0;

% 系统状态方程系数
af = [1 -2 1] .* (dT^-2) + [0 kk/m 0] + [1 0 -1] ./ 2 / dT * c / m;
% 系统输入系数
bf = 1 / m;

% 使用状态方程模拟船只位置
xtrue = filter(bf, af, wk);
xpred = rand(M,1);
% 使用一阶差分方程模拟船只速度
xptrue = filter([1 -1] / dT, 1, xtrue);

% 计算测量值
sigma = 0.3;     % 测量噪声标准差
a = 0.2;         % 海底坡度常数
b = 0;
z = sin(xtrue) + a * xtrue + b * xtrue.^2 + randn(M, 1) * sqrt(sigma); % 测量值，包含海底和噪声

% 初始化粒子
xk = rand(N, 1) * (xmax - xmin) + xmin;  % 初始化粒子位置
xpk = zeros(N, 1);                        % 初始化粒子速度
pik = repmat(1/N, N, 1);                   % 初始化粒子权重

% 绘图初始化
figure(3)
set(3, 'doublebuffer', 'on', 'position', [239 291 681 343])
t = linspace(xmin, xmax, N)';
h1 = line(t, sin(t) + a*t + b*t.^2, 'marker', 'none'); % 海底曲线
h3 = line(xtrue(1), 0, 'marker', 'o', 'linestyle', 'none'); % 真实船只位置
boat.x = [-2.4 -2 -1 0 1 2 2 -2.4]'; % 船只形状
boat.y = [1 0 -0.2 -0.3 -0.3 -0.3 0.8 1]' + 5; % 船只形状
hboat = line(boat.x + xtrue(1), boat.y, 'color', 'r'); % 船只位置
hecho = line(xtrue(1) * [1 1]', [z(1) boat.y(4)]', 'color', 'r', 'linestyle', '--'); % 测量值线
hwater = line([xmin xmax], boat.y(2) * [1 1], 'color', 'b', 'linestyle', '--'); % 水面线
h_hist = line(1, 1, 'color', 'k'); % 权重分布直方图
xlabel('Position x (m)')
xlim([xmin, xmax])
ylim([-3 ceil(max(boat.y)) + 2])
htext(1) = text(xmin + 0.5, boat.y(2), 'Sea surface', 'verticalalignment', 'bottom'); % 文本：海面
htext(2) = text(boat.x(7) + xtrue(1), boat.y(2), 'S/Y OptFilt', 'verticalalignment', 'bottom', 'horizontalalignment', 'right'); % 文本：船只名称
htext(3) = text(xmin + 0.5, 0, 'p(x)', 'verticalalignment', 'bottom'); % 文本：概率密度函数
htext(4) = text(xmin + 0.5, sin(xmin + 0.5) + a * (xmin + 0.5) + b * (xmin + 0.5).^2, 'Sea bottom', 'verticalalignment', 'bottom', 'color', 'b'); % 文本：海底


resample = 0.5;
Neff = zeros(M, 1);

sigma_sqrt_2_pi = sigma * sqrt(2 * pi);
two_sigma_square = 2 * sigma^2;

% 主循环
for k = 1:M
    % 时间更新
    wk = randn(N, 1) * F0;
    xk = xk + xpk * dT;
    xpk = xpk + (wk - xk * kk - xpk * (c - dT * kk)) / m * dT;
    
    % 测量更新
    pik = pik .* exp(-(sin(xk) + a * xk + b * xk.^2 - z(k)).^2 / two_sigma_square) / sigma_sqrt_2_pi;
    pik = pik / sum(pik);
    
    % 需要时进行重采样
    if resample > 0
        Neff(k) = 1 / sum(pik.^2);
        if Neff(k) < (resample * N)
            Inew = rsmp(pik, N);
            xk = xk(Inew);
            xpk = xpk(Inew);
            pik = repmat(1/N, N, 1);
        end
    end
    
    % [max_pik, max_index] = max(pik);
    % xpred(k) = xk(max_index);
    xpred(k) = dot(pik,xk);
    
    % 更新绘图
    set(h3, 'xdata', xtrue(k), 'ydata', z(k))
    set(hboat, 'xdata', boat.x + xtrue(k));
    set(hecho, 'xdata', xtrue(k) * [1 1]', 'ydata', [z(k) boat.y(4)]');
    
    [plotx, ploty] = histweight(xk, pik, 200, [xmin, xmax]);
    set(h_hist, 'xdata', plotx, 'ydata', ploty);
    set(htext(2), 'position', [boat.x(7) + xtrue(k), boat.y(2)]);
    
    drawnow;
end

% 绘制粒子效率曲线
figure(4)
plot(Neff/N)
xlabel('time step')
ylabel('Particle efficiency')
title('Efficient number of particles')

figure()
plot(xtrue, 'b'); hold on;
plot(z, 'r'); hold on;
plot(xpred, 'g');
legend('True Position', 'Measurements', 'Particle Filtered Estimate');
xlabel('Time Step');
ylabel('Position');
title('Particle Filter Prediction of Boat Position');

