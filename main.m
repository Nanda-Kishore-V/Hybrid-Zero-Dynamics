precomputation()

import autogen_func.swing_model
import autogen_func.impact_model
import autogen_func.hzd_LgLfy
import autogen_func.hzd_Lf2y

q0_minus = [pi/8, -pi/8, 0.4965, 1.5535, -2.0370, -0.4325]';
q0_plus = impact_model(q0_minus(1), q0_minus(2), q0_minus(3), q0_minus(4), q0_minus(5), q0_minus(6));

load('+data/hd.mat');

H = [1 0 -1; 0 1 -1; 1 0 0];
b = [pi; pi; 0];

q0_minus = [H*q0_minus(1:3) + b; H*q0_minus(4:6)];
q0_plus = [H*q0_plus(1:3) + b; H*q0_plus(4:6)];

load('+data/V_body.mat', 'V_body');
load('+data/D_body.mat', 'D_body');

gamma_0 = D_body(end, :);

qb1 = sym('qb1', 'real');
qb2 = sym('qb2', 'real');
qN = sym('qN', 'real');

qb = [qb1; qb2; qN];

dqb1 = sym('dqb1' ,'real');
dqb2 = sym('dqb2' ,'real');
dqN = sym('dqN' ,'real');

K = 0.5*[dqb1, dqb2, dqN]*D_body*[dqb1; dqb2; dqN];

epsilon1 = q0_plus(3);
epsilon2 = double(subs(jacobian(K, dqN), [qb1, qb2, qN, dqb1, dqb2, dqN], q0_plus'));

theta = qN;
hd1 = poly2sym(hd(1, :), qN);
hd2 = poly2sym(hd(2, :), qN);
h1 = qb1 - hd1;
h2 = qb2 - hd2;
h = [h1; h2];
k1 = jacobian(theta, qb)*inv([jacobian(h, qb); gamma_0])*[0; 0; 1];
k1 = subs(k1, [qb1, qb2], [hd1, hd2]);
k2 = -jacobian(V_body, qN);
k2 = subs(k2, [qb1, qb2], [hd1, hd2]);

V0 = -int(k2/k1, qN, [-pi/8, qN]);
V0 = matlabFunction(V0);
V0_minus = V0(pi/8);
arg_V0_max = fminbnd(@(x) -V0(x), -pi/8, pi/8);
V0_max = V0(arg_V0_max);

k1_func = simplify(k1);
k2_func = simplify(k2);
k1 = matlabFunction(k1);
k2 = matlabFunction(k2);
% figure(1)
% fplot(k1, [-pi/8, pi/8]);
% figure(2)
% fplot(k2, [-pi/8, pi/8]);
% return

% lambda_dq = [jacobian(h, qb); gamma_0]\[0; 0; 1]; %\
% lambda_dq = double(subs(lambda_dq, [qb1, qb2, qN], q0_minus(1:3)'));
% 
% q1 = sym('q1', 'real');
% q2 = sym('q2', 'real');
% q3 = sym('q3', 'real');
% sym('p0h', 'real');
% sym('p0v', 'real');
% q = H\([qb1; qb2; qN] - b);
% qbe = [qb1; qb2; qN; p0h; p0v];
% 
% inv(jacobian(H*[q1; q2; q3]+b, qb))'*De*inv(jacobian(q, qs))
% De = subs(De, [q1 q2 q3], q');
% P2 = subs(P2, [q1 q2 q3], q');
% E2 = simplify(jacobian(P2, qbe));
% del_F2 = -(E2*(De\E2'))\E2*[eye(3);zeros(2,3)];
% del_qedot = De\E2'*del_F2 + [eye(3);zeros(2,3)];
% 
% del_qdot = subs(del_qdot, [q1 q2 q3], q');
% del_qdot = double(subs(del_qdot, [qb1, qb2, qN], q0_minus(1:3)'));
% del_zero = double(subs(gamma_0, [qb1, qb2, qN], q0_plus(1:3)'))*del_qdot*lambda_dq;
del_zero = double(subs(jacobian(K, dqN), [qb1, qb2, qN, dqb1, dqb2, dqN], q0_plus'))/double(subs(jacobian(K, dqN), [qb1, qb2, qN, dqb1, dqb2, dqN], q0_minus'));

epsilon2_minus = epsilon2/del_zero;
zeta2_minus = 0.5*epsilon2_minus^2;
zeta2_plus = del_zero^2*zeta2_minus;

epsilon2_func = int(k2_func/k1_func, qN, [-pi/8, qN]);
epsilon2_func = matlabFunction(epsilon2_func);
int_function = @(x) 1/(k1(x)*sqrt(epsilon2^2 + 2*epsilon2_func(x)));
T = integral(int_function, -0.3927, 0.3927, 'ArrayValued', true);

odefunc = @(~, y) [k1(y(1))*y(2); k2(y(1))];

qbs = [];
times = [];
us = [];
total_impacts = 3;
start_time = 0;
for i = 1:total_impacts    
    options = odeset('Events', @eventsfn);
    [t, y, te, ye, ie] = ode45(@(t, y) odefunc(t, y), [start_time:0.01:start_time+10,Inf], [epsilon1; epsilon2], options);
    
    % impact
    epsilon1 = -pi/8;
    epsilon2 = del_zero*y(end, 2);
    
    start_time = t(end);
    
    for j = y'
        qb = [polyval(hd(1,:), j(1)); ...
              polyval(hd(2,:), j(1)); ...
              j(1)];

        dqb = [polyval(polyder(hd(1,:)), j(1)); ...
               polyval(polyder(hd(2,:)), j(1)); ...
               k1(j(1))*j(2)];

        q = H\(qb - b);
        dq = H\dqb;

        u = -inv(hzd_LgLfy(q(1), q(2), q(3), dq(1), dq(2), dq(3))) ...
            *hzd_Lf2y(q(1), q(2), q(3), dq(1), dq(2), dq(3));

        qbs = [qbs; q'];
        us = [us; u'];
    end
    times = [times; t];
end
disp('The torque squared cost is:');
disp(trapz(times', sum(us.^2, 2)')/total_impacts);
figure(1);
plot(times, qbs);
figure(2);
plot(times, us);

function [value, isterminal, direction] = eventsfn(~, y)
    value = y(1) >= pi/8;
    isterminal = 1;
    direction  = 0;
end