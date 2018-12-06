import autogen_func.swing_model
import autogen_func.impact_model

q0_minus = [0.3927, -0.3926, 0.4965, 1.5570, -2.0085, -0.4396]';
q0_plus = impact_model(q0_minus(1), q0_minus(2), q0_minus(3), q0_minus(4), q0_minus(5), q0_minus(6));

load('+data/hd.mat');

H = [1 0 -1; 0 1 -1; 1 0 0];
b = [pi; pi; 0];

q0_minus = [H*q0_minus(1:3) + b; H*q0_minus(4:6)];
q0_plus = [H*q0_plus(1:3) + b; H*q0_plus(4:6)];


[D, V] = change_coords(H, b);

gamma_0 = D(end, :);

qb1 = sym('qb1', 'real');
qb2 = sym('qb2', 'real');
qN = sym('qN', 'real');

q = [qb1; qb2; qN];

dqb1 = sym('dqb1' ,'real');
dqb2 = sym('dqb2' ,'real');
dqN = sym('dqN' ,'real');

K = 0.5*[dqb1, dqb2, dqN]*D*[dqb1; dqb2; dqN];

epsilon1 = q0_plus(3);
epsilon2 = double(subs(jacobian(K, dqN), [qb1, qb2, qN, dqb1, dqb2, dqN], q0_plus'))

theta = qN;
hd1 = poly2sym(hd(1, :), qN);
hd2 = poly2sym(hd(2, :), qN);
h1 = qb1 - hd1;
h2 = qb2 - hd2;
h = [h1; h2];
k1 = jacobian(theta, q)*inv([jacobian(h, q); gamma_0])*[0; 0; 1];
k1 = subs(k1, [qb1, qb2], [hd1, hd2]);
k2 = -jacobian(V, qN);
k2 = subs(k2, [qb1, qb2], [hd1, hd2]);

k1 = matlabFunction(k1);
k2 = matlabFunction(k2);
% figure(1)
% fplot(k1, [-pi/8, pi/8]);
% figure(2)
% fplot(k2, [-pi/8, pi/8]);
% return

% lambda_dq = [jacobian(h, q); gamma_0]\[0; 0; 1];
% lambda_dq = double(subs(lambda_dq, [qb1, qb2, qN], q0_minus(1:3)'));
% 
% [~, del_qdot] = compute_impact_model();
% q1 = sym('q1', 'real');
% q2 = sym('q2', 'real');
% q3 = sym('q3', 'real');
% q = H\([qb1; qb2; qN] - b);
% del_qdot = subs(del_qdot, [q1 q2 q3], q');
% del_qdot = double(subs(del_qdot, [qb1, qb2, qN], q0_minus(1:3)'));
% del_zero = double(subs(gamma_0, [qb1, qb2, qN], q0_plus(1:3)'))*del_qdot*lambda_dq;
del_zero = double(subs(jacobian(K, dqN), [qb1, qb2, qN, dqb1, dqb2, dqN], q0_plus'))/double(subs(jacobian(K, dqN), [qb1, qb2, qN, dqb1, dqb2, dqN], q0_minus'));

e1 = [];
e2 = [];
t = [];
dt = 0.01;
for i = 0:dt:10
    if epsilon1 >= pi/8
        epsilon1 = -pi/8;
        epsilon2 = del_zero*epsilon2;
    end
    epsilon1 = epsilon1 + k1(epsilon1)*epsilon2*dt;
    epsilon1 = wrapToPi(epsilon1);
    epsilon2 = epsilon2 + k2(epsilon1)*dt;
    t = [t, i];
    e1 = [e1, epsilon1];
    e2 = [e2, epsilon2];
end
figure(1);
plot(t, e1);
figure(2);
plot(t, e2);
e1(end)
e2(end)