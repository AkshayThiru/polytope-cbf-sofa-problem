function [u, h, dhdt_e, log] = control(i, t, s, env, param)
    % Environment.
    A = env.A;
    b = env.b;
    C = env.C;
    d = env.d;
    
    % Hyper-parameters.
    a1    = param.a1;
    P     = param.P;
    a2    = param.a2;
    eps_2 = param.eps_2;
    v_M   = param.v_M;
    w_M   = param.w_M;
    head  = param.head;
    sf    = param.sf;
    
    % 
    p = s(1:2);
    R = [cos(s(3)) -sin(s(3)); sin(s(3)) cos(s(3))];
    dR = [-sin(s(3)) -cos(s(3)); cos(s(3)) -sin(s(3))];
    
    mo_i = [1 5 9 13];
    mr_i = [1 5 9];
    
    % Evaluate h between all pairs of robots.
    lmb = zeros(8,6);
    Ndist = zeros(1,6);
    
    options_q = optimoptions('quadprog','Display','off','LinearSolver','sparse',...
        'MaxIterations',5000);
    
    log.time  = zeros(1, 7);
    log.duals = zeros(8, 6);
    
    for o = 1:3
        for r = 1:2
            k = (o-1)*2 + r;
            
            Ai = A(mo_i(o):mo_i(o+1)-1,:);
            bi = b(mo_i(o):mo_i(o+1)-1,1);
            Ci = C(mr_i(r):mr_i(r+1)-1,:)*R';
            di = d(mr_i(r):mr_i(r+1)-1,1) + Ci*p;
            
            Aeq = [Ci' Ai'];
            beq = zeros(2,1);
            Ain = -eye(8);
            bin = zeros(8,1);
            H = 1/4*[Ci*Ci' zeros(4,4);
                     zeros(4,4) Ai*Ai'];
            f = [di; bi];
            tic;
            [sol,fval] = quadprog(H,f,Ain,bin,Aeq,beq,[],[],[], options_q);
            log.time(2*(o-1)+r) = toc;
            lmb(:,k) = sol(1:8);
            Ndist(k) = fval;
        end
    end
    
    log.duals(:,:) = lmb;
    
    % Non-negativity constraints.
    temp = reshape(lmb,48,1);
    Z = eye(48);
    Z = -Z(temp < eps_2,:);
    
    Aeq = zeros(12,2+48+1);
    beq = zeros(12,1);
    Ain = zeros(1+6+size(Z,1)+4+1,2+48+1);
    bin = zeros(1+6+size(Z,1)+4+1,1);
    
    % CLF constraints
    Ain(1,:) = [(s(1:2)-sf(1:2))'*[cos(s(3)+head); sin(s(3)+head)] (s(3)-sf(3)) zeros(1,48) -1];
    bin(1,1) = -a1*((s-sf)'*(s-sf));
    
    Ain(end,:) = [zeros(1,2+48) -1];
    
    % Input bounds.
    Ain(1+6+size(Z,1)+1:1+6+size(Z,1)+4,1:2) = [eye(2); -eye(2)];
    bin(1+6+size(Z,1)+1:1+6+size(Z,1)+4,1) = [v_M; w_M; v_M; w_M];
    
    Ain(1+6+1:1+6+size(Z,1),:) = [zeros(size(Z,1),2) Z zeros(size(Z,1),1)];
    
    % Pairwise CBF constraints.
    for o = 1:3
        for r = 1:2
            k = 2*(o-1) + r;
            
            Ai = A(mo_i(o):mo_i(o+1)-1,:);
            bi = b(mo_i(o):mo_i(o+1)-1,1);
            Ci = C(mr_i(r):mr_i(r+1)-1,:);
            di = d(mr_i(r):mr_i(r+1)-1,1) + Ci*R'*p;
            
            l1 = lmb(1:4,k);
            l2 = lmb(5:8,k);
            
            Aeq(2*k-1:2*k,:) = [zeros(2,1) dR*Ci'*l1 zeros(2,8*(k-1)) R*Ci' Ai' zeros(2,8*(6-k)+1)];
            Ain(k+1,:) = [l1'*Ci*R'*[cos(s(3)+head); sin(s(3)+head)] l1'*Ci*dR'*p zeros(1,8*(k-1)) di' bi'+0.5*l2'*(Ai*Ai') zeros(1,8*(6-k)+1)];
            bin(k+1,1) = -2*a2*(Ndist(1,k) + 0.015^2);
        end
    end
    
    % Cost function.
    H = [eye(2,2) zeros(2,48+1); zeros(48,2+48+1); zeros(1,2+48) P] + 1e-2*eye(2+48+1);
    f = zeros(2+48+1,1);
    
    x0 = zeros(2+48+1,1);
    x0(end) = 10;
    
    % Solve problem.
    options.cl = [-inf*ones(1+6+size(Z,1)+4+1,1); zeros(12,1)];   % Lower bounds on the constraint functions.
    options.cu = [zeros(1+6+size(Z,1)+4+1,1); zeros(12,1)];   % Upper bounds on the constraint functions.
    funcs.constraints = @(x) [Ain; Aeq]*x - [bin; beq];
    funcs.jacobian = @(x) sparse([Ain; Aeq]);
    funcs.jacobianstructure = @() sparse([Ain; Aeq]);
    
    funcs.objective = @(x) 0.5*x'*H*x + f'*x;
    funcs.gradient = @(x) H*x + f;
    funcs.hessian = @(x, sigma, lambda) sparse(sigma*H); 
    funcs.hessianstructure  = @() sparse(H);
    
    options.ipopt.print_level = 0;
    
    tic;
    [sol, ~] = ipopt(x0, funcs, options);
    log.time(7) = toc;
    
    % Extract outputs.
    u = sol(1:2);
    h = -Ndist;
    dhdt_e = -Ain(2:7,:)*sol;
    dhdt_e = dhdt_e';
end

