function [] = record_data()
%     curr_x = [-pi/8; pi/8; pi/6; 0.94; -0.28; 2.04];
    curr_x = [pi/8, -pi/8, pi/6, 1.6, -1.6, 0];
    curr_x = impact_model(curr_x(1), curr_x(2), curr_x(3), curr_x(4), curr_x(5), curr_x(6));
    
    start_time = 0;
    overall_y = [];
    overall_t = [];
    overall_u = [];
    total_time = 10;
    for a = 1:total_time
        options = odeset('Events', @eventsfn);
        [t, y, te, ye, ie] = ode45(@(t, y) odefunc(t, y), [start_time, start_time+10], curr_x, options);
        if ie ~= 1
            disp("The robot fell.");
            return
        end
        for b = y'
            x_cell = num2cell(b);
            control = LgLfy(x_cell{:})\(psi_func(x_cell{:}) - Lf2y(x_cell{:}));
            overall_u = [overall_u; [control(1), control(2)]];
        end
        overall_y = [overall_y; y];
        overall_t = [overall_t; t];
        start_time = t(end);
        x_cell = num2cell(y(end, :));
        curr_x = impact_model(x_cell{:});
    end
    q1 = wrapToPi(y(:, 1) - y(:, 3) + pi);
    q2 = wrapToPi(y(:, 2) - y(:, 3) + pi);
    th = wrapToPi(y(:, 1));
    disp("The torque squared cost is:");
    disp(trapz(overall_t', sum(overall_u.^2, 2)')/total_time);
    theta = sym('theta', 'real');
    poly1 = polyfit(th, q1, 7);
    poly2 = polyfit(th, q2, 7);
    hd = [poly1; poly2];
    save('hd.mat', 'hd');
    figure(1);
    plot(th, polyval(hd(1, :), th));
    hold on;
    plot(th, q1);
    hold off;
    figure(2);
    plot(th, polyval(hd(2, :), th));
    hold on;
    plot(th, q2);
    hold off;
end

function [dx] = odefunc(~, x)
    x_cell = num2cell(x);
    u = LgLfy(x_cell{:})\(psi_func(x_cell{:}) - Lf2y(x_cell{:}));
    dx = swing_model(x_cell{:}, u(1), u(2));
end

function [value, isterminal, direction] = eventsfn(~, y)
    value = [y(1) > pi/8, abs(y(1)) > pi/2 || abs(y(2)) > pi/2];%sum(abs(y(1:3)) > pi/2)];
    isterminal = [1, 1];
    direction  = [0, 0];
end