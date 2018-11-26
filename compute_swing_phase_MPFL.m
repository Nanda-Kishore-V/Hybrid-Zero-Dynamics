function [] = compute_swing_phase_MPFL()
    % constants
    l = sym('0.5');
    r = sym('1.0');
    M_T = sym('10.0');
    M_H = sym('15.0');
    m = sym('5.0');
    g = sym('9.81');
    
    % variables
    q1 = sym('q1', 'real');
    q2 = sym('q2', 'real');
    q3 = sym('q3', 'real');
    
    dq1 = sym('dq1', 'real');
    dq2 = sym('dq2', 'real');
    dq3 = sym('dq3', 'real');
    
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

    % link CoMs
    G1 = P0 + r*e1/2;
    G2 = P1 + r*e2/2;
    
    % robot CoM
    m_tot = m+m+M_H+M_T;
    Pcm = (m*G1 + m*G2 + M_H*P1 + M_T*P3)/m_tot;
    
    q = [q1; q2; q3];
    qdot = [dq1; dq2; dq3];
    
    dPcm = jacobian(Pcm, q)*qdot;
    
    % Modified coordinates of the form
    % x_tilde = (qb; theta; dqb; sigma_N)
    qb1 = sym('qb1', 'real');
    qb2 = sym('qb2', 'real');
    theta = sym('theta', 'real');
    dqb1 = sym('dqb1', 'real');
    dqb2 = sym('dqb2', 'real');
    sigma_N = sym('sigma_N', 'real');
    
    dtheta = sym('dtheta', 'real');
    
    q_tilde = [qb1; qb2; theta];
    qb = [qb1; qb2];
    dqb = [dqb1; dqb2];
    
    x_tilde = [qb1; qb2; theta; dqb1; dqb2; sigma_N];
    
    % Modified controls
    v1 = sym('v1', 'real');
    v2 = sym('v2', 'real');
    v = [v1; v2];
    
    % Canonical change of coordinates
    % Note to self: This change of coordinate works because differentiation
    % of a constant is zero.
    H = [1 0 -1; 0 1 -1; 1 0 0];
    b = [pi;pi; 0];
    
    % swing phase model
    [D, C, G, B, V, K] = compute_swing_model();
    
    q = H\(q_tilde - b);
    V_tilde = subs(V, [q1 q2 q3], [q(1), q(2), q(3)]);
    V_tilde = simplify(V_tilde);
    D_tilde = inv(H)'*D*inv(H);
    D_tilde = subs(D_tilde, [q1 q2 q3], [q(1), q(2), q(3)]);
    D_tilde = simplify(D_tilde);
    
    J_norm = D_tilde(end, 1:end-1)/D_tilde(end, end);
    
    dx_tilde = [dqb; (sigma_N/D_tilde(end, end)) - J_norm*dqb; v; -diff(V_tilde, theta)];
    
    load('hd.mat');
    hd1 = poly2sym(hd(1,:), theta);
    hd2 = poly2sym(hd(2,:), theta);
    hd = [hd1; hd2];
    
    dhdth = jacobian(hd, theta); % dh_d / d theta
    LgLfh_tilde = eye(2) + dhdth*J_norm;
    LgLfh_tilde_inv = eye(2) + (dhdth*J_norm)/(1 + J_norm*dhdth);
    ddtheta = -diff(V_tilde, theta)/D_tilde(end, end) - (sigma_N/D_tilde(end, end)^2)*jacobian(D_tilde(end,end), qb)*dqb - jacobian(J_norm*dqb, qb)*dqb;
    dtheta = sigma_N/D_tilde(end, end) - J_norm*dqb;
    Lf2h_MPFL = - dhdth*ddtheta - jacobian(dhdth, theta)*dtheta^2;
    
%     matlabFunction(LgLfh_MPFL, 'File', 'LgLfh_MPFL', 'vars', [qb1, qb2, theta, dqb1, dqb2, sigma_N]);
%     matlabFunction(dx_tilde, 'File', 'swing_phase_MPFL', 'vars', [qb1, qb2, theta, dqb1, dqb2, sigma_N, v1, v2]);
    
    I = D_tilde(end, end) + D_tilde(end, 1:end-1)*dhdth;
    I = subs(I, [qb1 qb2], [hd1 hd2]);
    
    k1 = 1.0/I;
    q = H\([hd1; hd2; theta] - b);
    k2 = m_tot*g*subs(Pcm, [q1 q2 q3], q');
    k2 = k2(1);
    
    figure(1);
    fplot(k1, [-0.4, 0.4]);
    figure(2);
    fplot(k2, [-0.4, 0.4]);
    
    epsilon1 = -pi/8;
    epsilon2 = double(subs(D(end, 1:end), [q1 q2 q3], [-pi/8, pi/8, pi/6])*[0.9366; -0.2755; 2.0323]);
    
    dt = 0.01;
    e1 = [];
    t = [];
    for i = 0:0.01:1
        t = [t, i];
        e1 = [e1, epsilon1];
        epsilon1 = epsilon1 + double(subs(k1, theta, epsilon1))*epsilon2*dt;
        epsilon2 = epsilon2 + double(subs(k2, theta, epsilon1))*dt;
        if epsilon1 >= pi/8
            break
        end
    end
    double(subs(hd1, theta, epsilon1))
    double(subs(hd2, theta, epsilon1))
    epsilon1
    plot(t, e1);
    
%     dPcmv = subs(dPcm, [q1 q2 q3], q');
%     substitution_vars = H\[diff(hd1, theta)*dq1; diff(hd2, theta)*dq1; dq1];
%     dPcmv = subs(dPcmv, [dq1 dq2 dq3], substitution_vars');
% %     dPcmv(2)
% %     I*dq1
% %     simplify(dPcmv(2)/(I*dq1))
% 
%     % https://arxiv.org/pdf/1706.01127.pdf eqn 64
%     d = P2(2);
%     dzero = 1 + m_tot*d*diff(dPcmv(2), theta)*k1;
%     vpa(dzero, 2)
    
end