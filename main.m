%% Set path.
path(path, genpath('./'));

%% Environment and parameter setup.
% Environment polytope: Ai.z <= bi, i = {1,2,3}.
rect  = [eye(2); -eye(2)];
env.A = [rect; rect; rect];
env.b = [[0; 7; 2; -0.5]; [8; 7; 2; -6]; [8; 5; -1; -0.5]];
% Sofa polytope: Ci.z <= di, i = {1,2}.
env.C = [rect; rect];
env.d = [[0.1; 0; 0; 0.9+0.1]; [0.9+0.1; 0; 0; 0.1]];

% Hyper-parameters.
param.a1    = 0.1; % CLF constant [1/s].
param.a2    = 1; % CBF constant [1/s].
param.eps_2 = 1e-5; % U.S.C. parameter [].
param.v_M   = 0.3; % Velocity bounds [m/s].
param.w_M   = 0.2; % Angular velocity bounds [rad/s].
param.P     = 10; % CLF slack variable weight.

% Initial and final states.
param.head = pi/4; % Initial heading angle [rad].
param.sf   = [7; 5.95; -pi/4]; % Final state of sofa.
s0 = [0.05; 1.5; param.head]; % Initial state: position [m, m], orientation [rad].

% Display options.
display_solve_stats = true;
plot_stat_hist      = false;
display_animation   = false;
save_animation      = false;
plot_snapshot       = true;
save_snapshot       = false;
print_figures       = false;

%% Run simulation.
T  = 30; % Simulation time [s].
dt = 0.05; % Sample time [s].
tspan = [0 T];

% Logging data structures.
t      = (0:dt:T+dt)'; % Time.
s      = zeros(size(t,1),3); % State.
h      = zeros(size(t,1),6); % Square of minimum distance [m^2].
dhdt_e = zeros(size(t,1),6); % Estimated derivative of h [m^2/s].
s(1,:) = s0';

log.time  = zeros(size(t,1)-1,7); % Computation time per loop.
log.duals = zeros(8,6,size(t,1)-1); % Dual variables.

% Sample-and-hold simulation.
display_text = 0;
for i = 1:length(t)-1
    % Compute control input.
    tic
    [u_t,h(i,:),dhdt_e(i,:), log_i] = control(i, dt*(i-1), s(i,:)', env, param);
    loop_time = toc;
    % Statistics.
    fprintf(repmat('\b', 1, display_text));
    display_text = fprintf("time: %2.2f s, loop time: %1.4f s, frequency: %3.1f Hz", ...
        t(i), loop_time, 1/loop_time);
    % Simulate system.
    [t_te, s_te] = ode45(@(t,s) dyn_u(t,s,u_t,param), [0 dt], s(i,:)');
    s(i+1,:) = s_te(end,:);
    % Logging.
    log.time(i,:)    = log_i.time;
    log.duals(:,:,i) = log_i.duals;
end
fprintf('\n');

%% Statistics, plotting, and animations.
% Solve time statistics.
pctl = @(v,p) interp1(linspace(0.5/length(v), 1-0.5/length(v), length(v))', sort(v), p*0.01, 'spline');

if display_solve_stats
    fprintf('min. dist QPs: mean = %1.4f ms, std = %1.4f ms, p50 = %1.4f ms, p99 = %1.4f ms\n', ...
        mean(mean(log.time(:,1:6))), std(log.time(1:size(log.time,1)*6)), ...
        pctl(log.time(1:size(log.time,1)*6), 50), pctl(log.time(1:size(log.time,1)*6), 99));
    
    fprintf('control QP   : mean = %1.4f ms, std = %1.4f ms, p50 = %1.4f ms, p99 = %1.4f ms\n', ...
        mean(log.time(:,7)), std(log.time(:,7)), pctl(log.time(:,7), 50), pctl(log.time(:,7), 99));   
    
    fprintf('Total        : mean = %1.4f ms, std = %1.4f ms, p50 = %1.4f ms, p99 = %1.4f ms\n', ...
        mean(sum(log.time,2)), std(sum(log.time,2)), pctl(sum(log.time,2), 50), pctl(sum(log.time,2), 99));
end

if plot_stat_hist
    figure;
    histogram(time(:,1:6));
    figure;
    histogram(time(:,7));
end

% Animation.
if display_animation
	animate(t, s, 250, env);
end

if save_animation
	movie(t, s, 15, 25, 0);
end

% Snapshot image.
if plot_snapshot
	snapshot(t, s, [4.5, 7, 10, 12, 14, 16, 17.8, 20, 24, 30], env, save_snapshot);
end




% Plot trajectory of the vertex of the sofa.
figure;
hold on;
plot(s(:,1),s(:,2),'-b');
grid on
hold off;

% Set colours.
color_band = [[0, 0.4470, 0.7410]
    [0.8500, 0.3250, 0.0980];
    [0.9290, 0.6940, 0.1250];
    [0.4940, 0.1840, 0.5560];
    [0.4660, 0.6740, 0.1880];
    [0.6350, 0.0780, 0.1840];
    ];

fig_L = 600;
fig_H = 200;

% Plot minimum distance function.
figure('Renderer', 'painters', 'Position', [0 0 4/3*fig_L fig_H]);
til = tiledlayout(1, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
for arm = 1:2
    nexttile;
    hold on
    plot(t,h(:,1+(arm-1)), '-', 'Color', color_band(1,:), 'LineWidth', 1.5);
    plot(t,h(:,3+(arm-1)), '-', 'Color', color_band(2,:), 'LineWidth', 1.5);
    plot(t,h(:,5+(arm-1)), '-', 'Color', color_band(3,:), 'LineWidth', 1.5);
    ax = gca;
    ax.XTick = 0:10:30;
    ax.YTick = [10^-4, 10^-3, 10^-2, 10^-1, 1, 10];
    xlim([0 30]);
    if arm == 1
        ylim([10^-4 100]);
    else
        set(gca, 'YTickLabel', []);
    end
    box on
    set(gca,'LineWidth', 1, 'FontSize', 10);
    set(gca, 'YScale', 'log');
    h_l=get(gca,'Children');
    h_legend = legend(h_l,...
        {strcat('$h^{',num2str(arm),'3}(t)$'), strcat('$h^{',num2str(arm),'2}(t)$'), strcat('$h^{',num2str(arm),'1}(t)$')});
    h_legend.FontSize = 15;
    h_legend.ItemTokenSize = [15, 18];
    h_legend.NumColumns = 1;
    set(h_legend, 'Interpreter','latex');
    hold off
end
xlabel(til, 'Time (s)','interpreter','latex','FontSize', 15);
ylabel(til, '(m$^2$)','interpreter','latex','FontSize', 15);
if print_figures
	print(gcf,strcat('./img/sofa-h.png'), '-dpng', '-r1000');
	print(gcf,strcat('./img/sofa-h.eps'), '-depsc');
end

% Plot CBF enforcement.
arm = 1;
figure('Renderer', 'painters', 'Position', [0 0 4/3*fig_L fig_H]);
til = tiledlayout(1, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
for wall = 1:2:3
    nexttile;
    k = 2*(wall-1) + arm;
    hold on
    plot(t,-2*param.a2*h(:,k), 'LineWidth', 1.5);
    dhdt = gradient(h(:,k),t);
    plot(t(1:end-2),dhdt(1:end-2), 'LineWidth', 1.5);
    plot(t,dhdt_e(:,k), 'LineWidth', 1.5);
    ax = gca;
    ax.XTick = 0:10:30;
    if k == 1
        ax.YTick = -6:2:2;
        ylim([-6 2]);
    elseif k == 5
        ax.YTick = -0.15:0.05:0.1;
        ylim([-0.15 0.1]);
    end
    xlim([0 30]);
    box on
    set(gca,'LineWidth', 1, 'FontSize', 10);
    h_l=get(gca,'Children');
    armStr = num2str(arm);
    wallStr = num2str(wall);
    h_legend = legend(h_l,...
        {strcat('$\hat{\dot{h}}^{',armStr,wallStr,'}(t)$'), strcat('$\dot{h}^{',armStr,wallStr,'}(t)$'), strcat('$\bar{\alpha}^{',armStr,wallStr,'}(t)$')});
    h_legend.FontSize = 15;
    h_legend.ItemTokenSize = [15, 18];
    h_legend.NumColumns = 1;
    set(h_legend, 'Interpreter','latex');
    hold off
end
xlabel(til, 'Time (s)','interpreter','latex','FontSize', 15);
ylabel(til, '(m$^2$/s)','interpreter','latex','FontSize', 15);
if print_figures
    print(gcf,strcat('./img/sofa-dhdt.png'), '-dpng', '-r1000');
    print(gcf,strcat('./img/sofa-dhdt.eps'), '-depsc');
end

% Plot dual variables.
wall = 3;
figure('Renderer', 'painters', 'Position', [0 0 4/3*fig_L fig_H]);
til = tiledlayout(1, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
for arm = 1:2
    nexttile;
    k = 2*(wall-1) + arm;
    hold on
    idx = [1, 4, 6, 7];
    for i = 1:4
        plot(t(1:end-1),squeeze(log.duals(idx(i),k,:)), '-', 'Color', color_band(i,:), 'LineWidth', 1.5);
    end
    ax = gca;
    ax.XTick = 0:10:30;
    xlim([0 30]);
    ylim([0 0.7]);
    if arm == 2
        set(gca, 'YTickLabel', []);
    end
    box on
    set(gca,'LineWidth', 1, 'FontSize', 10);
    h_l=get(gca,'Children');
    armStr = num2str(arm);
    h_legend_str = {};
    for i = 1:4
        if i > 2
            h_legend_str{end+1} = strcat('$\lambda_',num2str(mod(idx(i)-1,4)+1),'^{A^',armStr,'}(t)$');
        else
            h_legend_str{end+1} = strcat('$\lambda_',num2str(mod(idx(i)-1,4)+1),'^{W^3}(t)$');
        end
    end
    h_legend = legend(h_l,h_legend_str);
    h_legend.FontSize = 15;
    h_legend.ItemTokenSize = [15, 18];
    h_legend.NumColumns = 1;
    set(h_legend, 'Interpreter','latex');
    % % grid on
    hold off
end
xlabel(til, 'Time (s)','interpreter','latex','FontSize', 15);
ylabel(til, 'Duals','interpreter','latex','FontSize', 15);
if print_figures
    print(gcf,strcat('./img/dual-soln.png'), '-dpng', '-r1000');
    print(gcf,strcat('./img/dual-soln.eps'), '-depsc');
end
