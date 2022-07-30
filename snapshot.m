function [] = snapshot(t, s, ts, env, save_fig)
    A = env.A;
    b = env.b;
    C = env.C;
    d = env.d;
    
    b = [[0; 7; 2; -2.0];
     [8; 7; 2; -6];
     [8; 5; -1; -2.0]];
    
    sa = spline(t,s',ts)';
    
    figure('Renderer', 'painters', 'Position', [0 0 600 300]);
    set(gca,'LooseInset',get(gca,'TightInset'));
    hold on
    
    % Plot environment.
    O1 = Polyhedron('A', A(1:4,:), 'b', b(1:4,1));
    O2 = Polyhedron('A', A(5:8,:), 'b', b(5:8,1));
    O3 = Polyhedron('A', A(9:12,:), 'b', b(9:12,1));
    O1.plot('color', 'red', 'linewidth', 1, 'linestyle', '-', 'wire', false, 'alpha', 1);
    O2.plot('color', 'red', 'linewidth', 1, 'linestyle', '-', 'wire', false, 'alpha', 1);
    O3.plot('color', 'red', 'linewidth', 1, 'linestyle', '-', 'wire', false, 'alpha', 1);
    
    text(-1,4.0,'$W_1$','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex');
    text(3,6.5,'$W_2$','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex');
    text(4.5,3.5,'$W_3$','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex');
    % Plot sofa.
    for i = 1:length(ts)
        R = [cos(sa(i,3)) -sin(sa(i,3)); sin(sa(i,3)) cos(sa(i,3))];
        P2 = Polyhedron('V',sa(i,1:2) + [0 0; [1 0]*R'; [0 -0.1]*R'; [1 -0.1]*R']);
        P1 = Polyhedron('V',sa(i,1:2) + [0 0; [0.1 0]*R'; [0 -1]*R'; [0.1 -1]*R']);

        P1.plot('color', 'b', 'linewidth', 1, 'linestyle', '-', 'wire', false, 'alpha', 1-i/12);
        P2.plot('color', 'g', 'linewidth', 1, 'linestyle', '-', 'wire', false, 'alpha', 1-i/12);
    end
    xlim([-2.0, 8.0]);
    ylim([2.0, 7.0]);
    set(gca,'xticklabel',[]);
    set(gca,'yticklabel',[]);
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    title('Snapshots of the moving sofa problem', 'interpreter', 'latex', 'FontSize', 15);
    axis('equal');
    box off
    grid off
    if save_fig
        print(gcf,'./img/non-conv-snapshot.png', '-dpng', '-r1000');
        print(gcf,'./img/non-conv-snapshot.eps', '-depsc', '-r600');
    end
end

