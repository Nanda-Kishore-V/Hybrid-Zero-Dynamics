function [D, C, G, B, V_tot, K_tot] = compute_swing_model()
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

    ddq1 = sym('ddq1', 'real');
    ddq2 = sym('ddq2', 'real');
    ddq3 = sym('ddq3', 'real');

    u1 = sym('u1', 'real');
    u2 = sym('u2', 'real');

    qs = [q1; q2; q3];
    qsdot = [dq1; dq2; dq3];
    qsddot = [ddq1; ddq2; ddq3];
    u = [u1; u2];

    load('+data/potential_energy.mat');
    load('+data/kinetic_energy.mat');
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

    xdot = [qsdot; D\(-C*qsdot - G + B*u)];
    matlabFunction(xdot, 'File', '+autogen_func/swing_model', 'vars', [q1 q2 q3 dq1 dq2 dq3 u1 u2]);
    save('+data/D.mat', 'D');
    save('+data/C.mat', 'C');
    save('+data/G.mat', 'G');
    save('+data/B.mat', 'B');    
end