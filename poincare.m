function [] = poincare()

    i = @(a) [pi/8, -pi/8, pi/6, a, -a, 0];
    
    dt1 = 1:0.1:6;
    
    rho = repelem(Inf, size(dt1, 2));
    for index = 1:size(dt1, 2)
        x_minus = i(dt1(index));
        xm_cell = num2cell(x_minus);
        x_plus = impact_model(xm_cell{:});

        options = odeset('Events', @eventsfn);
        [t, y, te, ye, ie] = ode45(@(t, y) odefunc(t, y), [0, 100], x_plus, options);
        if ie == 1 && prod(y(end, 1:3) - [pi/8, -pi/8, pi/6] < 0.01)
            rho(index) = ye(1, 4);
        end
    end
    plot(dt1, rho, 'LineWidth', 3);
    hold on;
    plot(dt1, dt1, 'LineWidth', 1);
end

function [dx] = odefunc(~, x)
    u = simple_controller(x(1), x(2), x(3), x(4), x(5), x(6));
    dx = swing_model(x(1), x(2), x(3), x(4), x(5), x(6), u(1), u(2));
end

function [value, isterminal, direction] = eventsfn(~, y)
    value = [y(1) > pi/8, sum(abs(y(1:3)) > pi/2)];
    isterminal = [1, 1];
    direction  = [0, 0];
end