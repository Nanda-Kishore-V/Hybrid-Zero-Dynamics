function [] = demo_script(controller, savevideo)
    import autogen_func.impact_model
    import autogen_func.swing_model
    if controller == 3
        import autogen_func.controller_ctl
    elseif controller == 2
        import autogen_func.fc_LgLfy
        import autogen_func.fc_Lf2y
        import autogen_func.fc_psi_func
    else
        import autogen_func.hzd_LgLfy
        import autogen_func.hzd_Lf2y
        import autogen_func.hzd_Lfy
        import autogen_func.hzd_y
    end
    
    swing_phase_model = @(q1,q2,q3,dq1,dq2,dq3,u1,u2) swing_model(q1,q2,q3,dq1,dq2,dq3,u1,u2); 
    simple_controller = @(q1,q2,q3,dq1,dq2,dq3) controller_ctl(q1,q2,q3,dq1,dq2,dq3);
    LgLfy = @(q1,q2,q3,dq1,dq2,dq3) fc_LgLfy(q1,q2,q3,dq1,dq2,dq3);
    Lf2y = @(q1,q2,q3,dq1,dq2,dq3) fc_Lf2y(q1,q2,q3,dq1,dq2,dq3);
    psi_func = @(q1,q2,q3,dq1,dq2,dq3) fc_psi_func(q1,q2,q3,dq1,dq2,dq3);
    LgLfy_hzd = @(q1,q2,q3,dq1,dq2,dq3) hzd_LgLfy(q1,q2,q3,dq1,dq2,dq3);
    Lf2y_hzd = @(q1,q2,q3,dq1,dq2,dq3) hzd_Lf2y(q1,q2,q3,dq1,dq2,dq3);
    Lfy_hzd = @(q1,q2,q3,dq1,dq2,dq3) hzd_Lfy(q1,q2,q3,dq1,dq2,dq3);
    y_hzd = @(q1,q2,q3,dq1,dq2,dq3) hzd_y(q1,q2,q3,dq1,dq2,dq3);
    
    curr_x = [pi/8, -pi/8, pi/6, 1.6, -1.6, 0];
    curr_x = impact_model(curr_x(1), curr_x(2), curr_x(3), curr_x(4), curr_x(5), curr_x(6));
    stance_leg_end = [0, 0];
    
    start_time = 0;
    overall_y = [];
    overall_t = [];
    overall_u = [];
    total_impacts = 10;
    
    if savevideo
        writerObj = VideoWriter('biped_walking.avi');
        writerObj.FrameRate = 30;
        open(writerObj);
    end
    
    epsilon = 0.1;
    Kd = eye(2);
    Kp = eye(2);
    for a = 1:total_impacts
        options = odeset('Events', @eventsfn);
        [t, y, te, ye, ie] = ode45(@(t, y) odefunc(t, y), [start_time:0.01:start_time+10, Inf], curr_x, options);
        if ie ~= 1
            disp('Debugging:'); disp(te); disp(ie);  disp(ye);
            simulate(y, stance_leg_end);
            disp("The robot fell in step " + num2str(a));
            return
        end
        u_timestep = [];
        for b = y'
            x_cell = num2cell(b);
            if controller == 1
                control = -LgLfy_hzd(x_cell{:})\(Lf2y_hzd(x_cell{:}) + (1/epsilon)*Kd*Lfy_hzd(x_cell{:}) + (1/epsilon^2)*Kp*y_hzd(x_cell{:}));
            elseif controller == 2
                control = LgLfy(x_cell{:})\(psi_func(x_cell{:}) - Lf2y(x_cell{:}));
            else
                control = simple_controller(x_cell{:});
            end
            u_timestep = [u_timestep; [control(1), control(2)]];
        end
        frames = simulate(y, stance_leg_end);
        
        if savevideo
            for frame = frames
                writeVideo(writerObj, frame);
            end
        end
        
        overall_y = [overall_y; y];
        overall_t = [overall_t; t];
        overall_u = [overall_u; u_timestep];
        start_time = t(end);
        x_cell = num2cell(y(end, :));
        curr_x = impact_model(x_cell{:});
        stance_leg_end = [stance_leg_end(1) + sin(curr_x(2)) - sin(curr_x(1)); 0];
        disp("The torque squared cost in this orbit is:");
        disp(trapz(t', sum(u_timestep.^2, 2)'));
    end
    if savevideo
        close(writerObj);
    end
    overall_y(end ,:)
    figure(2);
    plot(overall_t, overall_y(:,1:3));
    hold on;
    grid on;
    xlabel('$t$','Interpreter','latex')
    ylabel('$\theta_1, \theta_2, \theta_3$','Interpreter','latex')
    title('Configuration variables  vs time','Interpreter','latex')
    leg1 = legend('$\theta_1$', '$\theta_2$', '$\theta_3$');
    set(leg1,'Interpreter','latex');
    axis tight;
    hold off;
    figure(3);
    plot(overall_t, overall_u);
    hold on;
    grid on;
    xlabel('$t$','Interpreter','latex')
    ylabel('$u_1, u_2$','Interpreter','latex')
    title('Controls  vs time','Interpreter','latex')
    leg1 = legend('$u_1$','$u_2$');
    set(leg1, 'Interpreter', 'latex');
    axis tight;
    hold off;
    disp("The torque squared cost is:");
    disp(trapz(overall_t', sum(overall_u.^2, 2)')/total_impacts);

    function [dx] = odefunc(~, x)
        x(1:3) = wrapToPi(x(1:3));
        x_c = num2cell(x);
        if controller == 1
            u = -LgLfy_hzd(x_c{:})\(Lf2y_hzd(x_c{:}) + (1/epsilon)*Kd*Lfy_hzd(x_c{:}) + (1/epsilon^2)*Kp*y_hzd(x_c{:}));
        elseif controller == 2
            u = LgLfy(x_c{:})\(psi_func(x_c{:}) - Lf2y(x_c{:}));
        else
            u = simple_controller(x_c{:});
        end
        dx = swing_phase_model(x_c{:}, u(1), u(2));
    end

    function [value, isterminal, direction] = eventsfn(~, y)
        value = [y(1) >= pi/8, abs(y(1)) > pi/2 || abs(y(2)) > pi/2];%sum(abs(y(1:3)) > pi/2)];
        isterminal = [1, 1];
        direction  = [0, 0];
    end
end

function [frames] = simulate(y, stance_leg_end)
    iter = 1;
    frames = [];
%     for a = y'
%         qe = [a(1:3); stance_leg_end(1); stance_leg_end(2)];
%         frames(iter) = draw_robot(qe);
%         iter = iter + 1;
%         pause(0.0001);
%     end
end