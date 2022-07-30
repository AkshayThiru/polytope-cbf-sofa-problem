function [] = animate(t, s, N, env)
    figure;
    
    A = env.A;
    b = env.b;
    C = env.C;
    d = env.d;
    
    ta = linspace(0,t(end),N);
    sa = spline(t,s',ta)';
    
    for i = 1:length(ta)
        % Construct polytopes.
        R = [cos(sa(i,3)) -sin(sa(i,3)); sin(sa(i,3)) cos(sa(i,3))];
        P1 = Polyhedron('A', C(1:4,:)*R', 'b', d(1:4,1)+C(1:4,:)*R'*sa(i,1:2)');
        P2 = Polyhedron('A', C(5:8,:)*R', 'b', d(5:8,1)+C(5:8,:)*R'*sa(i,1:2)');
        O1 = Polyhedron('A', A(1:4,:), 'b', b(1:4,1));
        O2 = Polyhedron('A', A(5:8,:), 'b', b(5:8,1));
        O3 = Polyhedron('A', A(9:12,:), 'b', b(9:12,1));
        hold on
        P1.plot();
        P2.plot();
        O1.plot();
        O2.plot();
        O3.plot();
        xlim([-4, 6]);
        ylim([0, 10]);
        text(-1,8,['time = ',num2str(i/length(ta)*t(end))]);
        hold off
        pause(0.01);
        if i ~= length(ta)
            clf
        end
    end
end

