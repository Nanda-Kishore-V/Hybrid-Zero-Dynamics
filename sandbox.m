f = figure();
% curr_x = [-pi/8; pi/8; pi/6; 0.94; -0.28; 2.04];
curr_x = [pi/8, -pi/8, pi/6, 1.6, -1.6, 0];

stance_leg_end = [0; 0];

odefunc = @(t, x, u1, u2) swing_model(x(1), x(2), x(3), x(4), x(5), x(6), u1, u2);

x1 = [];    x2 = [];    x3 = [];
u1 = [];    u2 = [];

sim_time = 0:0.01:5;

J = 0;

figure(1);
title('Simulation of bipedal robot');

for t = sim_time
    x1 = [x1, curr_x(1)];
    x2 = [x2, curr_x(2)];
    x3 = [x3, curr_x(3)];
    u1 = [u1, 0];
    u2 = [u2, 0];
    
    x_cell = num2cell(curr_x);
    if abs(curr_x(1) - pi/8) < 0.1
        disp(J);
        curr_x = impact_model(x_cell{:});
        stance_leg_end = [stance_leg_end(1) + sin(curr_x(2)) - sin(curr_x(1)); 0];
        qe = [curr_x(1:3); stance_leg_end];
        draw_robot(qe);
        J = 0;
        continue
    end
%     u = simple_controller(x_cell{:});
    u = LgLfy(x_cell{:})\(psi_func(x_cell{:}) - Lf2y(x_cell{:}));
    J = J + sum(u.^2)*0.01;
    u1(end) = u(1);
    u2(end) = u(2);
    u_cell = num2cell(u);
    dx = swing_model(x_cell{:}, u_cell{:});
    [~, x] = ode45(@(t, y) odefunc(t, y, u_cell{:}), [t, t+0.01], curr_x);
    curr_x = x(end,:)';
    
    qe = [curr_x(1:3); stance_leg_end];
    draw_robot(qe);
    pause(0.0001);
end 

figure(2);
x = [x1; x2; x3];
plot(sim_time, x);
leg1 = legend('$\theta_{1}$','$\theta_{2}$','$\theta_{3}$');
set(leg1,'Interpreter','latex');
title('States');

figure(3);
u = [u1; u2];
plot(sim_time, u);
leg1 = legend('$u_{1}$','$u_{2}$');
set(leg1,'Interpreter','latex');
title('Control inputs');