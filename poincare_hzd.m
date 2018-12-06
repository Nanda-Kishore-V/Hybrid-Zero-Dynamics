function [] = poincare_hzd()
    import autogen_func.swing_model
    import autogen_func.impact_model
    
    import autogen_func.hzd_LgLfy
    import autogen_func.hzd_Lfy
    import autogen_func.hzd_Lf2y
    import autogen_func.hzd_y
    
    swing_phase_model = @(q1,q2,q3,dq1,dq2,dq3,u1,u2) swing_model(q1,q2,q3,dq1,dq2,dq3,u1,u2); 
    LgLfy_hzd = @(q1,q2,q3,dq1,dq2,dq3) hzd_LgLfy(q1,q2,q3,dq1,dq2,dq3);
    Lf2y_hzd = @(q1,q2,q3,dq1,dq2,dq3) hzd_Lf2y(q1,q2,q3,dq1,dq2,dq3);
    Lfy_hzd = @(q1,q2,q3,dq1,dq2,dq3) hzd_Lfy(q1,q2,q3,dq1,dq2,dq3);
    y_hzd = @(q1,q2,q3,dq1,dq2,dq3) hzd_y(q1,q2,q3,dq1,dq2,dq3);
    
    load('+data/hd.mat', 'hd');
    
    i = @(a) [pi/8; -pi/8; pi/6; a; -a; 0];
    
    dq1 = 1:0.1:10;
    rho = repelem(Inf, size(dq1, 2));
    
    for index = 1:size(dq1, 2)
        x_minus = i(dq1(index));
        xm_cell = num2cell(x_minus);
        x_plus = impact_model(xm_cell{:});
%         xb_minus = [H*x_minus(1:3) + b; H*x_minus(4:6)];
%         xbm_cell = num2cell(xb_minus);
%         xb_plus = [H*x_plus(1:3) + b; H*x_plus(4:6)];
%         xbp_cell = num2cell(xb_plus);
%         z_minus = [x_minus(1); momentum_conj(xbm_cell{:})];
%         z_plus = [x_plus(1); momentum_conj(xbp_cell{:})];
        
        options = odeset('Events', @eventsfn);
        [t, y, te, ye, ie] = ode45(@(t, y) odefunc(t, y), [0, Inf], x_plus, options);
        dtheta = ye(1,4);
        for idx = 1:15
            if ie ~= 1
                break
            end
            xm_cell = num2cell(y(end,:));
            x_plus = impact_model(xm_cell{:});
            [~, y, ~, ye, ie] = ode45(@(t, y) odefunc(t, y), [0, 3], x_plus, options);
        end
        if ie == 1
            rho(index) = dtheta;
        end
    end
    plot(dq1, rho, 'LineWidth', 3);
    hold on;
    grid on;
    plot(dq1, dq1, 'LineWidth', 1);
    xlabel('$\dot{\theta_1}^{-}$','Interpreter','latex')
    ylabel('$\rho(\dot{\theta_1}^{-})$','Interpreter','latex')
    title('Poincare map')
    
    function [dx] = odefunc(~, x)
        epsilon = 0.1;
        Kd = eye(2);
        Kp = eye(2);
        x_c = num2cell(x);
        u = -LgLfy_hzd(x_c{:})\(Lf2y_hzd(x_c{:}) + (1/epsilon)*Kd*Lfy_hzd(x_c{:}) + (1/epsilon^2)*Kp*y_hzd(x_c{:}));
        dx = swing_phase_model(x_c{:}, u(1), u(2));
    end

    function [value, isterminal, direction] = eventsfn(~, y)
        value = [y(1) >= pi/8, abs(y(1)) > pi/2 || abs(y(2)) > pi/2];
        isterminal = [1, 1];
        direction  = [0, 0];
    end
end