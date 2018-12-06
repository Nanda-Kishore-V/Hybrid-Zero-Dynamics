function [] = generate_feedback_controller()
% variables
q1 = sym('q1', 'real');
q2 = sym('q2', 'real');
q3 = sym('q3', 'real');

dq1 = sym('dq1', 'real');
dq2 = sym('dq2', 'real');
dq3 = sym('dq3', 'real');

qs = [q1; q2; q3];
qsdot = [dq1; dq2; dq3];

load('+data/D.mat', 'D');
load('+data/C.mat', 'C');
load('+data/G.mat', 'G');
load('+data/B.mat', 'B');

% x_dot = f(x) + g(x)u
f = [qsdot; D\(-C*qsdot - G)];
g = [zeros(3,2); D\B];

hd1 = 0.512 + 0.073*q1 + 0.035*q1^2 - 0.819*q1^3;
hd2 = -q1 + (-2.27 + 3.26*q1 + 3.11*q1^2 + 1.89*q1^3)*(q1 + pi/8)*(q1 - pi/8);

% desired output
y = simplify([q3 - hd1; q2 - hd2]);

Lfy = simplify(jacobian(y, [qs; qsdot])*f);
LgLfy = simplify(jacobian(Lfy, [qs; qsdot])*g);
Lf2y = simplify(jacobian(Lfy, [qs; qsdot])*f);

psi = simplify(psi_f(y, Lfy));
matlabFunction(LgLfy, 'File', '+autogen_func/fc_LgLfy', 'vars', [q1 q2 q3 dq1 dq2 dq3]);
matlabFunction(psi, 'File', '+autogen_func/fc_psi_func', 'vars', [q1 q2 q3 dq1 dq2 dq3]);
matlabFunction(Lf2y, 'File', '+autogen_func/fc_Lf2y', 'vars', [q1 q2 q3 dq1 dq2 dq3]);

    function [v] = psi_f(y, dy)
        epsilon = 0.1;
        alpha = 0.99;    

        phi_a = @(x1, x2) x1 + (sign(x2)*abs(x2)^(2-alpha))/(2-alpha);
        psi_a = @(x1, x2) -1*sign(x2)*abs(x2)^alpha - sign(phi_a(x1, x2))*abs(phi_a(x1, x2))^(alpha/(2-alpha));

        v = [psi_a(y(1), epsilon*dy(1))/epsilon^2;
             psi_a(y(2), epsilon*dy(2))/epsilon^2];
    end

end