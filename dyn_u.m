function [ds] = dyn_u(~, s, u, param)
    % Unicycle dynamics
    ds = [u(1)*cos(s(3)+param.head);
          u(1)*sin(s(3)+param.head);
          u(2)];
end

