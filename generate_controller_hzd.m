function [] = generate_controller_hzd()
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

load('+data/hd.mat', 'hd');
hd1 = poly2sym(hd(1,:), q1);
hd2 = poly2sym(hd(2,:), q1);

% desired output
y = simplify([q1 - q3 + pi - hd1; q2 - q3 + pi - hd2]);

Lfy = simplify(jacobian(y, [qs; qsdot])*f);
LgLfy = simplify(jacobian(Lfy, [qs; qsdot])*g);
Lf2y = simplify(jacobian(Lfy, [qs; qsdot])*f);

matlabFunction(Lfy, 'File', '+autogen_func/hzd_Lfy', 'vars', [q1 q2 q3 dq1 dq2 dq3]);
matlabFunction(LgLfy, 'File', '+autogen_func/hzd_LgLfy', 'vars', [q1 q2 q3 dq1 dq2 dq3]);
matlabFunction(Lf2y, 'File', '+autogen_func/hzd_Lf2y', 'vars', [q1 q2 q3 dq1 dq2 dq3]);
matlabFunction(y, 'File', '+autogen_func/hzd_y', 'vars', [q1 q2 q3 dq1 dq2 dq3]);
end
