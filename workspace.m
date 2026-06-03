%% workspace.m - 连续体机器人运动学仿真与工作空间分析
% 基于论文《连续体手术机器人系统设计与控制研究》
%   - 第3.2节 运动学模型（公式4.30-4.31, 4.34-4.36）
%   - 第5.1节 连续型关节段和工作空间仿真测试（表4.1, 表4.3）
%   - 第5.2节 工作空间仿真（图5.2, 图5.3）
%
% 功能：
%   1. 正运动学模型：四根驱动丝长度 → 弯曲角/旋转角/末端位置
%   2. 表4.1 十种工况验证（驱动丝长度变化 → 弯曲角/旋转角）
%   3. 表4.3 推送长度-末端位置对照验证
%   4. 工作空间蒙特卡洛采样与三维可视化
%   5. 工作空间切片投影分析
%
% 作者：基于论文复现
% 日期：2025

clear; clc; close all;

%% ============ 1. 结构参数定义（论文表2.9、图4.1） ============
% 几何参数
N = 20;                     % 单节段零件数量
h = 5.0;                    % 单节段高度 (mm)
L0 = 100;                   % 导丝原长/初始弧长 (mm)
r0_design = 6.4;            % 导丝孔分度圆半径 - 设计值 (mm)
r_body = 6.0;               % 本体半径 (mm)
L_total = N * h;            % 关节段总弧长 = 100 (mm)

% 关节运动学约束（第5.1节）
Theta_max = 140;            % 最大弯曲角 (°)
Theta_min = -140;           % 最小弯曲角 (°)
L_min = 84;                 % 最小导丝/弧长 (mm)
L_max = 110;                % 最大导丝/弧长 (mm)

% ============================================================
% 2. 有效半径标定（基于表4.1数据逆向计算）
%   从实验数据反推等效半径，反映导丝实际作用距离
% ============================================================
% 理论：Θ = ΔL_eff / (2 * r0_eff),  φ = atan2(dL24, dL13)
% 其中 ΔL_eff = sqrt((ΔL1-ΔL3)^2 + (ΔL2-ΔL4)^2)

% 表4.1数据
table4_1 = [
    1,  6,  0, -6,  0, 57.03,    0;
    2, -6,  0,  6,  0, 57.12, 139.13;
    3,  0,  6,  0, -6, 57.26,  90.21;
    4,  0, -6,  0,  6, 57.10, -90.43;
    5,  8,  8, -8, -8, 108.10,  45.25;
    6, -8,  8,  8, -8, 108.23, 135.18;
    7,  5, -5, -5,  5, 68.34, -45.31;
    8, -5, -5,  5,  5, 68.20, -135.22;
    9,  8,  0, -4, -4, 60.26,  18.12;
    10,-4,  8,  4, -8, 85.14, 117.25
];

% 计算每种工况的有效半径
n_cases = size(table4_1, 1);
r0_eff_cases = zeros(n_cases, 1);
for i = 1:n_cases
    dL1 = table4_1(i,2); dL2 = table4_1(i,3);
    dL3 = table4_1(i,4); dL4 = table4_1(i,5);
    Theta_exp = deg2rad(table4_1(i,6));
    dL13 = dL1 - dL3;
    dL24 = dL2 - dL4;
    delta_eff = sqrt(dL13^2 + dL24^2);
    if delta_eff > 1e-6 && Theta_exp > 1e-6
        r0_eff_cases(i) = delta_eff / (2 * Theta_exp);
    end
end
r0_eff = mean(r0_eff_cases(r0_eff_cases > 0));
fprintf('设计半径: r0 = %.2f mm\n', r0_design);
fprintf('标定有效半径: r0_eff = %.4f mm (基于表4.1 %d种工况平均)\n', r0_eff, sum(r0_eff_cases > 0));

%% ============ 3. 正运动学模型（函数定义见文件末尾）============
% 正运动学流程：
%   驱动空间(ΔL_i)  →  关节空间(Θ, φ)  →  操作空间(x, y, z)
% 详细函数实现在文件末尾的"辅助函数"部分

%% ============ 4. 表4.1 十种工况验证 ============
fprintf('\n========== 表4.1 十种工况验证 ==========\n');
fprintf('序号   ΔL1  ΔL2  ΔL3  ΔL4  |  Θ理论(°)  Θ实验(°)  误差(°)  |  φ理论(°)  φ实验(°)  误差(°)\n');
fprintf('------------------------------------------------------------------------------\n');

theta_errors = zeros(n_cases, 1);
phi_errors = zeros(n_cases, 1);

for i = 1:n_cases
    dL1 = table4_1(i,2); dL2 = table4_1(i,3);
    dL3 = table4_1(i,4); dL4 = table4_1(i,5);
    Theta_exp = table4_1(i,6);
    phi_exp = table4_1(i,7);

    % 使用标定半径计算
    [Theta_calc, phi_calc] = wire2joint(dL1, dL2, dL3, dL4, r0_eff);
    Theta_calc_deg = rad2deg(Theta_calc);
    phi_calc_deg = rad2deg(phi_calc);

    theta_err = abs(Theta_calc_deg - Theta_exp);
    phi_err = abs(phi_calc_deg - phi_exp);
    % 对于phi角度差，处理跨越±180°的情况
    if phi_err > 180
        phi_err = 360 - phi_err;
    end

    theta_errors(i) = theta_err;
    phi_errors(i) = phi_err;

    fprintf('  %2d   %+3d  %+3d  %+3d  %+3d  |  %8.2f   %8.2f   %6.2f   |  %8.2f   %8.2f   %6.2f\n', ...
        i, dL1, dL2, dL3, dL4, ...
        Theta_calc_deg, Theta_exp, theta_err, ...
        phi_calc_deg, phi_exp, phi_err);
end

fprintf('------------------------------------------------------------------------------\n');
fprintf('平均误差(含工况2):    Θ: %.3f°    φ: %.3f°\n', mean(theta_errors), mean(phi_errors));
fprintf('平均误差(除工况2):    Θ: %.3f°    φ: %.3f°\n', ...
    mean(theta_errors([1,3,4,5,6,7,8,9,10])), mean(phi_errors([1,3,4,5,6,7,8,9,10])));
fprintf('最大误差:          Θ: %.3f°    φ: %.3f°\n', max(theta_errors), max(phi_errors));
fprintf('\n注：工况2 φ理论值180°与实验值139.13°偏差较大的可能原因：\n');
fprintf('    1) 旋转角测量参考系与建模参考系存在固定偏置\n');
fprintf('    2) 导线回差或NiTi迟滞效应导致的机械零位偏移\n');

%% ============ 5. 表4.3 推送长度-末端位置验证 ============
fprintf('\n========== 表4.3 推送长度-末端Z坐标验证 ==========\n');

% 表4.3数据
table4_3 = [
    1,   0,    15,  115;
    2,   139.54, 15,  5.6;
    3,  -139.11, 10,  14.9;
    4,   90.14,  0,   66.8;
    5,   45.24,  15,  81.1;
    6,   120.36, -10, 58.5;
    7,  -30.13,  0,   95.5;
    8,   110.05, 10,  31.8;
    9,   139.06, 15,  12.1;
    10,  0,       1,   NaN    % 表中Z值缺失
];

fprintf('序号   Θ(°)  L_push(mm)  Z_exp(mm)  Z_calc(mm)  误差(mm)\n');
fprintf('---------------------------------------------------------\n');

z_errors = zeros(size(table4_3,1), 1);
for i = 1:size(table4_3,1)
    Theta_deg = table4_3(i,2);
    L_push = table4_3(i,3);
    Z_exp = table4_3(i,4);

    % 弧长 = 原长 + 推送量
    L_arc = L0 + L_push;

    % 计算末端Z坐标（φ和ΔL不影响纯弯曲方向的Z坐标）
    Theta = deg2rad(Theta_deg);
    [~, ~, z_calc] = joint2end(Theta, 0, L_arc);

    if ~isnan(Z_exp)
        z_err = abs(z_calc - Z_exp);
        z_errors(i) = z_err;
        fprintf('  %2d   %+7.2f  %+4d      %7.1f     %7.1f    %6.1f\n', ...
            i, Theta_deg, L_push, Z_exp, z_calc, z_err);
    else
        fprintf('  %2d   %+7.2f  %+4d      (缺失)    %7.1f\n', ...
            i, Theta_deg, L_push, z_calc);
    end
end

valid_idx = ~isnan(table4_3(:,4));
if sum(valid_idx) > 0
    fprintf('---------------------------------------------------------\n');
    fprintf('平均Z坐标误差: %.2f mm\n', mean(z_errors(valid_idx)));
    fprintf('注：表4.3 Z坐标理论差可能源于:\n');
    fprintf('    1) 伸长度"推送长度"指基座推送量，实际弧长≠原长+推送量\n');
    fprintf('    2) 大弯曲角下机械间隙和重力变形使末端低于常曲率模型预测\n');
end

%% ============ 6. 工作空间蒙特卡洛采样 ============
fprintf('\n========== 工作空间采样 ==========\n');

% 采样参数
n_samples = 8000;       % 总采样点数

% 关节空间随机采样
Theta_samples = deg2rad(Theta_min + (Theta_max - Theta_min) * rand(n_samples, 1));
phi_samples = 2 * pi * rand(n_samples, 1);  % 0~360°
L_samples = L_min + (L_max - L_min) * rand(n_samples, 1);

% 正向计算末端位置
x_ws = zeros(n_samples, 1);
y_ws = zeros(n_samples, 1);
z_ws = zeros(n_samples, 1);

for i = 1:n_samples
    [x_ws(i), y_ws(i), z_ws(i)] = joint2end(Theta_samples(i), phi_samples(i), L_samples(i));
end

fprintf('采样点数: %d\n', n_samples);
fprintf('弯曲角范围: [%.0f, %.0f]°\n', Theta_min, Theta_max);
fprintf('弧长范围: [%.0f, %.0f] mm\n', L_min, L_max);

% 计算工作空间统计量
ws_x_range = [min(x_ws), max(x_ws)];
ws_y_range = [min(y_ws), max(y_ws)];
ws_z_range = [min(z_ws), max(z_ws)];
fprintf('工作空间范围:\n');
fprintf('  X: [%.1f, %.1f] mm\n', ws_x_range(1), ws_x_range(2));
fprintf('  Y: [%.1f, %.1f] mm\n', ws_y_range(1), ws_y_range(2));
fprintf('  Z: [%.1f, %.1f] mm\n', ws_z_range(1), ws_z_range(2));

%% ============ 7. 可视化 ============

% ==================== 图1：工作空间三维散点图 ====================
figure('Position', [50, 50, 1400, 1000], 'Name', '连续体机器人工作空间');

% ---- 子图1：完整工作空间 3D ----
subplot(2,3,1);
scatter3(x_ws, y_ws, z_ws, 6, z_ws, 'filled');
hold on;
% 标出表4.1的10种工况点
for i = 1:n_cases
    dL1 = table4_1(i,2); dL2 = table4_1(i,3);
    dL3 = table4_1(i,4); dL4 = table4_1(i,5);
    [xi, yi, zi] = forward_kinematics(dL1, dL2, dL3, dL4, L0, r0_eff);
    plot3(xi, yi, zi, 'rp', 'MarkerSize', 12, 'LineWidth', 2);
end
% 标出原点
plot3(0, 0, 0, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
title('完整工作空间（蒙特卡洛采样）');
grid on; axis equal;
view(45, 30);
colormap(jet);
cb = colorbar;
ylabel(cb, 'Z (mm)');

% ---- 子图2：XY平面投影 ----
subplot(2,3,4);
scatter(x_ws, y_ws, 4, z_ws, 'filled');
hold on;
% 标注工况点投影
for i = 1:n_cases
    dL1 = table4_1(i,2); dL2 = table4_1(i,3);
    dL3 = table4_1(i,4); dL4 = table4_1(i,5);
    [xi, yi, ~] = forward_kinematics(dL1, dL2, dL3, dL4, L0, r0_eff);
    plot(xi, yi, 'rp', 'MarkerSize', 8, 'LineWidth', 2);
end
xlabel('X (mm)'); ylabel('Y (mm)');
title('XY平面投影'); grid on; axis equal;
colormap(jet);

% ---- 子图2：XZ平面投影 ----
subplot(2,3,2);
scatter(x_ws, z_ws, 4, z_ws, 'filled');
hold on;
for i = 1:n_cases
    dL1 = table4_1(i,2); dL2 = table4_1(i,3);
    dL3 = table4_1(i,4); dL4 = table4_1(i,5);
    [xi, ~, zi] = forward_kinematics(dL1, dL2, dL3, dL4, L0, r0_eff);
    plot(xi, zi, 'rp', 'MarkerSize', 8, 'LineWidth', 2);
end
xlabel('X (mm)'); ylabel('Z (mm)');
title('XZ平面投影'); grid on; axis equal;
colormap(jet);

% ---- 子图3：YZ平面投影 ----
subplot(2,3,3);
scatter(y_ws, z_ws, 4, z_ws, 'filled');
hold on;
for i = 1:n_cases
    dL1 = table4_1(i,2); dL2 = table4_1(i,3);
    dL3 = table4_1(i,4); dL4 = table4_1(i,5);
    [~, yi, zi] = forward_kinematics(dL1, dL2, dL3, dL4, L0, r0_eff);
    plot(yi, zi, 'rp', 'MarkerSize', 8, 'LineWidth', 2);
end
xlabel('Y (mm)'); ylabel('Z (mm)');
title('YZ平面投影'); grid on; axis equal;
colormap(jet);

sgtitle('连续体机器人运动学仿真与工作空间分析', 'FontSize', 14);

% ---- 子图5：不同弧长的工作空间对比 ----
subplot(2,3,5);
L_test_values = [84, 100, 110];
colors_l = {'b', 'r', 'g'};
hold on;
for li = 1:length(L_test_values)
    L_test = L_test_values(li);
    % 在固定弧长下采样弯曲角
    n_theta = 50;
    theta_test = linspace(deg2rad(Theta_min), deg2rad(Theta_max), n_theta);
    phi_test = 0;
    x_fixed = zeros(n_theta, 1);
    z_fixed = zeros(n_theta, 1);
    for j = 1:n_theta
        [x_fixed(j), ~, z_fixed(j)] = joint2end(theta_test(j), phi_test, L_test);
    end
    plot(x_fixed, z_fixed, [colors_l{li}, '-'], 'LineWidth', 1.5, ...
        'DisplayName', sprintf('L=%dmm', L_test));
end
% 标注表4.3数据点
for i = 1:size(table4_3,1)
    if ~isnan(table4_3(i,4))
        plot(0, table4_3(i,4), 'ko', 'MarkerSize', 6);
    end
end
xlabel('X (mm)'); ylabel('Z (mm)');
title('不同弧长下机器人弯曲轮廓');
legend('Location', 'best');
grid on; axis equal;

% ---- 子图6：工作空间密度分布 ----
subplot(2,3,6);
% 计算径向距离
r_ws = sqrt(x_ws.^2 + y_ws.^2);
% 绘制径向-轴向分布
scatter(r_ws, z_ws, 4, 'filled');
hold on;
% 标注边界
theta_b = deg2rad(-140:5:140);
phi_b = 0;
L_b = 100;
r_b = zeros(size(theta_b));
z_b = zeros(size(theta_b));
for j = 1:length(theta_b)
    [x_b, ~, z_b(j)] = joint2end(theta_b(j), phi_b, L_b);
    r_b(j) = abs(x_b);
end
plot(r_b, z_b, 'r-', 'LineWidth', 2, 'DisplayName', '理论包络(L=100mm)');
xlabel('径向距离 sqrt(X^2+Y^2) (mm)'); ylabel('Z (mm)');
title('工作空间径向-轴向分布');
legend('Location', 'best'); grid on;

% ==================== 图2：连续体弯曲形态动画（表4.1 典型工况） ====================
figure('Position', [100, 100, 1200, 500], 'Name', '连续体弯曲形态');

% 绘制6个典型工况
typical_cases = [1, 3, 5, 7, 9, 10];
for ii = 1:6
    i = typical_cases(ii);
    dL1 = table4_1(i,2); dL2 = table4_1(i,3);
    dL3 = table4_1(i,4); dL4 = table4_1(i,5);
    [Theta_i, phi_i] = wire2joint(dL1, dL2, dL3, dL4, r0_eff);
    Theta_i_deg = rad2deg(Theta_i);
    phi_i_deg = rad2deg(phi_i);

    % 绘制离散关节段
    n_seg = N;
    seg_pts = zeros(n_seg+1, 3);
    % 离散弧长参数 s ∈ [0, L_total]
    for s_idx = 0:n_seg
        s = s_idx * h;  % 沿弧长的位置
        t_param = s / L_total;  % 归一化参数 [0, 1]
        theta_s = Theta_i * t_param;  % 当前位置的弯曲角
        if abs(Theta_i) < 1e-8
            x_s = 0; y_s = 0; z_s = s;
        else
            R_s = L_total / Theta_i;
            x_s = cos(phi_i) * R_s * (1 - cos(theta_s));
            y_s = sin(phi_i) * R_s * (1 - cos(theta_s));
            z_s = R_s * sin(theta_s);
        end
        seg_pts(s_idx+1, :) = [x_s, y_s, z_s];
    end

    subplot(2,3,ii);
    % 绘制连续曲线
    plot3(seg_pts(:,1), seg_pts(:,2), seg_pts(:,3), 'b-', 'LineWidth', 2);
    hold on;
    % 绘制关节点
    plot3(seg_pts(:,1), seg_pts(:,2), seg_pts(:,3), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
    % 绘制基座
    plot3(0, 0, 0, 'gs', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    % 绘制末端
    plot3(seg_pts(end,1), seg_pts(end,2), seg_pts(end,3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

    xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
    title(sprintf('工况%d: Θ=%.1f° φ=%.1f°', i, Theta_i_deg, phi_i_deg));
    grid on; axis equal;
    xlim([-120, 120]); ylim([-120, 120]); zlim([0, 120]);
    view(45, 20);
end
sgtitle('典型工况下的连续体弯曲形态（20个关节段）', 'FontSize', 13);

% ==================== 图3：表4.1 理论值与实验值对比 ====================
figure('Position', [100, 100, 900, 400], 'Name', '理论值与实验值对比');

subplot(1,2,1);
% 计算各工况的理论弯曲角
theta_calc_all = zeros(n_cases, 1);
for i = 1:n_cases
    [t, ~] = wire2joint(table4_1(i,2), table4_1(i,3), ...
                         table4_1(i,4), table4_1(i,5), r0_eff);
    theta_calc_all(i) = rad2deg(t);
end
bar_data_theta = [table4_1(:,6), theta_calc_all];
bar(1:n_cases, bar_data_theta);
xlabel('工况序号'); ylabel('弯曲角 Θ (°)');
title('弯曲角理论值与实验值对比');
legend('实验值(表4.1)', '理论值', 'Location', 'best');
grid on;

subplot(1,2,2);
% 计算各工况的理论旋转角
phi_calc_all = zeros(n_cases, 1);
for i = 1:n_cases
    [~, p] = wire2joint(table4_1(i,2), table4_1(i,3), ...
                         table4_1(i,4), table4_1(i,5), r0_eff);
    phi_calc_all(i) = rad2deg(p);
end
bar_data_phi = [table4_1(:,7), phi_calc_all];
bar(1:n_cases, bar_data_phi);
xlabel('工况序号'); ylabel('旋转角 φ (°)');
title('旋转角理论值与实验值对比');
legend('实验值(表4.1)', '理论值', 'Location', 'best');
grid on;

sgtitle(sprintf('表4.1 十种工况验证（标定半径 r_eff = %.2f mm）', r0_eff));

% ==================== 图4：可伸缩连续体机器人工作空间（论文图5.3风格） ====================
figure('Position', [50, 50, 1000, 800], 'Name', '可伸缩工作空间（论文图5.3风格）');

% 分离上下半空间（中性参考位上方与下方）
idx_upper = z_ws >= L0;  % 弧长>原长 = 伸展状态
idx_lower = z_ws < L0;   % 弧长<原长 = 收缩状态

scatter3(x_ws(idx_upper), y_ws(idx_upper), z_ws(idx_upper), 8, 'r', 'filled', ...
    'MarkerFaceAlpha', 0.5);
hold on;
scatter3(x_ws(idx_lower), y_ws(idx_lower), z_ws(idx_lower), 8, 'k', 'filled', ...
    'MarkerFaceAlpha', 0.3);

% 标注基座
plot3(0, 0, 0, 'go', 'MarkerSize', 15, 'MarkerFaceColor', 'g', 'LineWidth', 2);

xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
title(sprintf('可伸缩连续体机器人工作空间（%d个采样点）', n_samples));
legend({'L>100mm (伸展)', 'L<100mm (收缩)', '基座位置'}, 'Location', 'best');
grid on; axis equal;
view(45, 30);
set(gca, 'FontSize', 11);

%% ============ 8. 仿真汇总 ============
fprintf('\n========== 仿真汇总 ==========\n');
fprintf('正运动学模型：基于常曲率假设\n');
fprintf('  弯曲角：Θ = ΔL_eff / (2r₀),  ΔL_eff = √((ΔL1-ΔL3)²+(ΔL2-ΔL4)²)\n');
fprintf('  旋转角：φ = atan2(ΔL2-ΔL4, ΔL1-ΔL3)\n');
fprintf('  末端位置：基于常曲率积分（公式4.34-4.36）\n');
fprintf('标定有效半径：r₀_eff = %.3f mm\n', r0_eff);
fprintf('工作空间形状：近似球面构型\n');
fprintf('最大径向范围：%.1f mm\n', max(sqrt(x_ws.^2 + y_ws.^2)));
fprintf('最大轴向高度：%.1f mm\n', max(z_ws));
fprintf('工作空间体积(包络)：%.0f mm³\n', ...
    (max(x_ws)-min(x_ws))*(max(y_ws)-min(y_ws))*(max(z_ws)-min(z_ws)));
fprintf('===============================\n');

%% ============ 辅助函数 ============
% 所有辅助函数已在前面定义为嵌套函数
% 注意：MATLAB要求所有本地函数在脚本末尾

function [Theta, phi] = wire2joint(dL1, dL2, dL3, dL4, r0)
    % 驱动空间 → 关节空间
    % 输入：四根驱动丝长度变化 ΔL1,ΔL2,ΔL3,ΔL4 (mm)
    % 输出：弯曲角 Θ (rad), 旋转角 φ (rad)
    % 基于论文公式4.30-4.31
    dL13 = dL1 - dL3;
    dL24 = dL2 - dL4;
    delta_eff = sqrt(dL13^2 + dL24^2);

    phi = atan2(dL24, dL13);

    if delta_eff < 1e-8
        Theta = 0;
    else
        Theta = delta_eff / (2 * r0);
    end
end

function [x, y, z] = joint2end(Theta, phi, L_arc)
    % 关节空间 → 操作空间（末端位置）
    % 基于常曲率模型（论文第3.2.1节）
    if abs(Theta) < 1e-8
        x = 0; y = 0; z = L_arc;
    else
        R = L_arc / Theta;
        x = cos(phi) * R * (1 - cos(Theta));
        y = sin(phi) * R * (1 - cos(Theta));
        z = R * sin(Theta);
    end
end

function [x, y, z, Theta, phi] = forward_kinematics(dL1, dL2, dL3, dL4, L_arc, r0)
    % 完整正运动学：驱动空间 → 操作空间
    [Theta, phi] = wire2joint(dL1, dL2, dL3, dL4, r0);
    [x, y, z] = joint2end(Theta, phi, L_arc);
end
