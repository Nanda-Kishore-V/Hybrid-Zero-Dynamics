function [] = generate_controller()
% constants
l = sym('0.5');
r = sym('1.0');
M_T = sym('10.0');
M_H = sym('15.0');
m = sym('5.0');
g = sym('9.81');

% constants
% l = sym('l', 'real');
% r = sym('r', 'real');
% M_T = sym('M_T', 'real');
% M_H = sym('M_H', 'real');
% m = sym('m', 'real');
% g = sym('g', 'real');


% variables
q1 = sym('q1', 'real');
q2 = sym('q2', 'real');
q3 = sym('q3', 'real');

dq1 = sym('dq1', 'real');
dq2 = sym('dq2', 'real');
dq3 = sym('dq3', 'real');

ddq1 = sym('ddq1', 'real');
ddq2 = sym('ddq2', 'real');
ddq3 = sym('ddq3', 'real');

u1 = sym('u1', 'real');
u2 = sym('u2', 'real');

% unit vectors
i = sym([1; 0]);
j = sym([0; 1]);
e1 = cos(-q1)*j - sin(-q1)*i;
e2 = -cos(-q2)*j + sin(-q2)*i;
e3 = cos(-q3)*j - sin(-q3)*i;

% link ends
P0 = 0*i + 0*j;
P1 = P0 + r*e1;
P2 = P1 + r*e2;
P3 = P1 + l*e3;

% link center of masses
G1 = P0 + r*e1/2;
G2 = P1 + r*e2/2;

% potential energies
V1 = m*g*G1(2);
V2 = m*g*G2(2);
VH = M_H*g*P1(2);
VT = M_T*g*P3(2);

qs = [q1; q2; q3];
qsdot = [dq1; dq2; dq3];
qsddot = [ddq1; ddq2; ddq3];
u = [u1; u2];

% derivative function
derivative = @(in) jacobian(in, qs)*qsdot;

% derivative of positions
dP1 = derivative(P1);
dP2 = derivative(P2);
dP3 = derivative(P3);
dG1 = derivative(G1);
dG2 = derivative(G2);

% kinetic energies
K1 = 0.5*m*dot(dG1, dG1);% + 0.5*m*(l^2)*(dq1^2)/12;
K2 = 0.5*m*dot(dG2, dG2);
K3 = 0.5*M_H*dot(dP1, dP1);
K4 = 0.5*M_T*dot(dP3, dP3);

V_tot = V1 + V2 + VH + VT;
K_tot = K1 + K2 + K3 + K4;

L = K_tot - V_tot;

dL_dq = simplify(jacobian(L, qs));
dL_dqdot = simplify(jacobian(L, qsdot))';
eqn = simplify([jacobian(dL_dqdot, qs), jacobian(dL_dqdot, qsdot)]*[qsdot;qsddot]);
eqn = simplify(eqn - dL_dq');
[D, ~] = equationsToMatrix(eqn, [ddq1, ddq2, ddq3]);
G = simplify(jacobian(V_tot, qs))';
C = sym(zeros(size(D,1), size(D,2)));
for k = 1:size(D,1)
    for j = 1:size(D,2)
        for i = 1:size(D,1)
            C(k,j) = C(k, j) + 0.5*(diff(D(k,j),qs(i)) + diff(D(k,i), qs(j)) - diff(D(i,j), qs(k)))*qsdot(i);
        end
    end
end
B = jacobian([q3 - q1; q3 - q2], qs)';
eqn = simplify(eqn == B*u);

% x_dot = f(x) + g(x)u
f = [qsdot; D\(-C*qsdot - G)];
g = [zeros(3,2); D\B];

% desired output
y = [q3 - pi/6; q2 + q1];

Lfy = simplify(jacobian(y, [qs; qsdot])*f);
LgLfy = simplify(jacobian(Lfy, [qs; qsdot])*g);
Lf2y = simplify(jacobian(Lfy, [qs; qsdot])*f);

u = LgLfy\(psi(y, Lfy) - Lf2y);
matlabFunction(u, 'File', 'simple_controller', 'vars', [q1 q2 q3 dq1 dq2 dq3]);

function [v] = psi(y, dy)
    epsilon = 0.1;
    alpha = 0.99;    
    
    phi_a = @(x1, x2) x1 + (sign(x2)*abs(x2)^(2-alpha))/(2-alpha);
    psi_a = @(x1, x2) -1*sign(x2)*abs(x2)^alpha - sign(phi_a(x1, x2))*abs(phi_a(x1, x2))^(alpha/(2-alpha));
    
    v = [psi_a(y(1), epsilon*dy(1))/epsilon^2;
         psi_a(y(2), epsilon*dy(2))/epsilon^2];
end

end