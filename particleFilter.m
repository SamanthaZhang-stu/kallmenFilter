%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SC23219041 ���� 2024/4/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ��������ռ䡢����ں�ͼ�δ���
clear, clc, close all

% ���崬ֻλ�õ�ȡֵ��Χ
xmin = -10;
xmax = 10;
% ���崬ֻ�ٶȵ�ȡֵ��Χ
xpmin = -0.3;
xpmax = 0.3;
% ������
N = 10000;
% ϵͳ����
m = 1;          % ������kg��
kk = 1;         % ���ɳ�����N/m��
c = 0.5;        % ����ϵ����N/s��
F0 = 10;        % �Ŷ���������ֵ��N��
dT = 0.05;      % ʱ�䲽����s��
M = 400;        % ʱ�䲽��

% ������ʵ���Ŷ�������
wk = (rand(M, 1) - 0.5) * 2 * F0;

% ϵͳ״̬����ϵ��
af = [1 -2 1] .* (dT^-2) + [0 kk/m 0] + [1 0 -1] ./ 2 / dT * c / m;
% ϵͳ����ϵ��
bf = 1 / m;

% ʹ��״̬����ģ�⴬ֻλ��
xtrue = filter(bf, af, wk);
xpred = rand(M,1);
% ʹ��һ�ײ�ַ���ģ�⴬ֻ�ٶ�
xptrue = filter([1 -1] / dT, 1, xtrue);

% �������ֵ
sigma = 0.3;     % ����������׼��
a = 0.2;         % �����¶ȳ���
b = 0;
z = sin(xtrue) + a * xtrue + b * xtrue.^2 + randn(M, 1) * sqrt(sigma); % ����ֵ���������׺�����

% ��ʼ������
xk = rand(N, 1) * (xmax - xmin) + xmin;  % ��ʼ������λ��
xpk = zeros(N, 1);                        % ��ʼ�������ٶ�
pik = repmat(1/N, N, 1);                   % ��ʼ������Ȩ��

% ��ͼ��ʼ��
figure(3)
set(3, 'doublebuffer', 'on', 'position', [239 291 681 343])
t = linspace(xmin, xmax, N)';
h1 = line(t, sin(t) + a*t + b*t.^2, 'marker', 'none'); % ��������
h3 = line(xtrue(1), 0, 'marker', 'o', 'linestyle', 'none'); % ��ʵ��ֻλ��
boat.x = [-2.4 -2 -1 0 1 2 2 -2.4]'; % ��ֻ��״
boat.y = [1 0 -0.2 -0.3 -0.3 -0.3 0.8 1]' + 5; % ��ֻ��״
hboat = line(boat.x + xtrue(1), boat.y, 'color', 'r'); % ��ֻλ��
hecho = line(xtrue(1) * [1 1]', [z(1) boat.y(4)]', 'color', 'r', 'linestyle', '--'); % ����ֵ��
hwater = line([xmin xmax], boat.y(2) * [1 1], 'color', 'b', 'linestyle', '--'); % ˮ����
h_hist = line(1, 1, 'color', 'k'); % Ȩ�طֲ�ֱ��ͼ
xlabel('Position x (m)')
xlim([xmin, xmax])
ylim([-3 ceil(max(boat.y)) + 2])
htext(1) = text(xmin + 0.5, boat.y(2), 'Sea surface', 'verticalalignment', 'bottom'); % �ı�������
htext(2) = text(boat.x(7) + xtrue(1), boat.y(2), 'S/Y OptFilt', 'verticalalignment', 'bottom', 'horizontalalignment', 'right'); % �ı�����ֻ����
htext(3) = text(xmin + 0.5, 0, 'p(x)', 'verticalalignment', 'bottom'); % �ı��������ܶȺ���
htext(4) = text(xmin + 0.5, sin(xmin + 0.5) + a * (xmin + 0.5) + b * (xmin + 0.5).^2, 'Sea bottom', 'verticalalignment', 'bottom', 'color', 'b'); % �ı�������


resample = 0.5;
Neff = zeros(M, 1);

sigma_sqrt_2_pi = sigma * sqrt(2 * pi);
two_sigma_square = 2 * sigma^2;

% ��ѭ��
for k = 1:M
    % ʱ�����
    wk = randn(N, 1) * F0;
    xk = xk + xpk * dT;
    xpk = xpk + (wk - xk * kk - xpk * (c - dT * kk)) / m * dT;
    
    % ��������
    pik = pik .* exp(-(sin(xk) + a * xk + b * xk.^2 - z(k)).^2 / two_sigma_square) / sigma_sqrt_2_pi;
    pik = pik / sum(pik);
    
    % ��Ҫʱ�����ز���
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
    
    % ���»�ͼ
    set(h3, 'xdata', xtrue(k), 'ydata', z(k))
    set(hboat, 'xdata', boat.x + xtrue(k));
    set(hecho, 'xdata', xtrue(k) * [1 1]', 'ydata', [z(k) boat.y(4)]');
    
    [plotx, ploty] = histweight(xk, pik, 200, [xmin, xmax]);
    set(h_hist, 'xdata', plotx, 'ydata', ploty);
    set(htext(2), 'position', [boat.x(7) + xtrue(k), boat.y(2)]);
    
    drawnow;
end

% ��������Ч������
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

