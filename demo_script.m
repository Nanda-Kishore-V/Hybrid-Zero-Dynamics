function [] = demo_script(controller)
    global controller_type;
    controller_type = controller;
    curr_x = [pi/8, -pi/8, pi/6, 1.6, -1.6, 0];
    curr_x = impact_model(curr_x(1), curr_x(2), curr_x(3), curr_x(4), curr_x(5), curr_x(6));
    stance_leg_end = [0, 0];
    
    start_time = 0;
    overall_y = [];
    overall_t = [];
    overall_u = [];
    total_time = 2;
    for a = 1:total_time
        options = odeset('Events', @eventsfn);
        [t, y, te, ye, ie] = ode45(@(t, y) odefunc(t, y), [start_time, start_time+10], curr_x, options);
        if ie ~= 1
            y
            for a = y'
                qe = [a(1:3); stance_leg_end(1); stance_leg_end(2)];
                draw_robot(qe);
                pause(0.0001);
            end
            disp("The robot fell.");
            return
        end
        for b = y'
            x_cell = num2cell(b);
            if controller == 1
                control = -LgLfy_hybrid(x_cell{:})\Lf2y_hybrid(x_cell{:});
            elseif controller == 2
                control = LgLfy(x_cell{:})\(psi_func(x_cell{:}) - Lf2y(x_cell{:}));
            else
                control = simple_controller(x_cell{:});
            end
            overall_u = [overall_u; [control(1), control(2)]];
        end
        for a = y'
            qe = [a(1:3); stance_leg_end(1); stance_leg_end(2)];
            draw_robot(qe);
            pause(0.0001);
        end

        overall_y = [overall_y; y];
        overall_t = [overall_t; t];
        start_time = t(end);
        x_cell = num2cell(y(end, :));
        curr_x = impact_model(x_cell{:});
        stance_leg_end = [stance_leg_end(1) + sin(curr_x(2)) - sin(curr_x(1)); 0];
    end
%     q1 = wrapToPi(y(:, 1) - y(:, 3) + pi);
%     q2 = wrapToPi(y(:, 2) - y(:, 3) + pi);
%     th = wrapToPi(y(:, 1));
    figure(2);
    plot(overall_t, overall_y(:,1:3));
    figure(3);
    plot(overall_t, overall_u);
    disp("The torque squared cost is:");
    disp(trapz(overall_t', sum(overall_u.^2, 2)')/total_time);
end

function [dx] = odefunc(~, x)
    global controller_type;
    x(1:3) = wrapToPi(x(1:3));
    x_cell = num2cell(x);
    if controller_type == 1
        u = -LgLfy_hybrid(x_cell{:})\Lf2y_hybrid(x_cell{:});
%         u = -LgLfy_hybrid(x_cell{:})\Lf2y_hybrid(x_cell{:});
    elseif controller_type == 2
        u = LgLfy(x_cell{:})\(psi_func(x_cell{:}) - Lf2y(x_cell{:}));
    else
        u = simple_controller(x_cell{:});
    end
    dx = swing_model(x_cell{:}, u(1), u(2));
end

function [value, isterminal, direction] = eventsfn(~, y)
    value = [y(1) > pi/8, abs(y(1)) > pi/2 || abs(y(2)) > pi/2];%sum(abs(y(1:3)) > pi/2)];
    isterminal = [1, 1];
    direction  = [0, 0];
end