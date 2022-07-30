function [] = movie(t, s, L, Tf, Ti, env)
    A = env.A;
    b = env.b;
    C = env.C;
    d = env.d;
    
    fig1 = figure();
    set(fig1, 'Position', [0 0 1550 1080]);
    xlim([-2.5, 8.5]);
    ylim([0, 7.5]);
    set(fig1, 'Renderer', 'painters');
    hold on
    
    rate = 50;
    N = rate*L;
    
    ta = linspace(Ti,Tf,N);
    sa = spline(t,s',ta)';
    
    v = VideoWriter('./mov/simulation.avi','Motion JPEG AVI');
    v.FrameRate = rate;
    v.Quality = 100;
    
    open(v);
    
    for i = 1:length(ta)
        R = [cos(sa(i,3)) -sin(sa(i,3)); sin(sa(i,3)) cos(sa(i,3))];
        
        P2 = Polyhedron('V',sa(i,1:2) + [0 0; [1 0]*R'; [0 -0.1]*R'; [1 -0.1]*R']);
        P1 = Polyhedron('V',sa(i,1:2) + [0 0; [0.1 0]*R'; [0 -1]*R'; [0.1 -1]*R']);
        
        O1 = Polyhedron('A', A(1:4,:), 'b', b(1:4,1));
        O2 = Polyhedron('A', A(5:8,:), 'b', b(5:8,1));
        O3 = Polyhedron('A', A(9:12,:), 'b', b(9:12,1));
        
        hold on

        P1.plot('color', 'b', 'linewidth', 1, 'linestyle', '-', 'wire', false);
        P2.plot('color', 'g', 'linewidth', 1, 'linestyle', '-', 'wire', false);
        O1.plot();
        O2.plot();
        O3.plot();
        
        text(4.5,1.5,['time = $',num2str((i-1)/(length(ta)-1)*Tf, '%.2f'), '$ $s$'],'FontSize',20,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex');
        text(-1,3.25,'$W_1$','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex');
        text(3,6.5,'$W_2$','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex');
        text(4.5,3.25,'$W_3$','FontSize',20,'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','latex');
        
        set(gca,'xticklabel',[]);
        set(gca,'yticklabel',[]);
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        set(gca,'visible','off');
        box off
        grid off
        
        drawnow;
        M = getframe(gcf);
        writeVideo(v,M);
        if i ~= length(ta)
            clf
        end
    end
    
    close(v);
end
