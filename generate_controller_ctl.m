function [] = generate_controller_ctl()
% ctl - constant torso lean
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

% desired output
y = [q3 - pi/6; q2 + q1];

Lfy = simplify(jacobian(y, [qs; qsdot])*f);
LgLfy = simplify(jacobian(Lfy, [qs; qsdot])*g);
Lf2y = simplify(jacobian(Lfy, [qs; qsdot])*f);

u = LgLfy\(psi(y, Lfy) - Lf2y);
matlabFunction(u, 'File', '+autogen_func/controller_ctl', 'vars', [q1 q2 q3 dq1 dq2 dq3]);

function [v] = psi(y, dy)
    epsilon = 0.1;
    alpha = 0.99;    
    
    phi_a = @(x1, x2) x1 + (sign(x2)*abs(x2)^(2-alpha))/(2-alpha);
    psi_a = @(x1, x2) -1*sign(x2)*abs(x2)^alpha - sign(phi_a(x1, x2))*abs(phi_a(x1, x2))^(alpha/(2-alpha));
    
    v = [psi_a(y(1), epsilon*dy(1))/epsilon^2;
         psi_a(y(2), epsilon*dy(2))/epsilon^2];
end

end