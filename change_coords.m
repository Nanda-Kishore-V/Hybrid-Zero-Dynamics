function [D_bar, V_bar] = change_coords(H, b)
    % variables
    q1 = sym('q1', 'real');
    q2 = sym('q2', 'real');
    q3 = sym('q3', 'real');

    dq1 = sym('dq1', 'real');
    dq2 = sym('dq2', 'real');
    dq3 = sym('dq3', 'real');
    
    qs = [q1; q2; q3];
    
    [D, C, G, B, V_tot, K_tot] = compute_swing_model();
    q_bar = H*qs + b;
    D_bar = inv(jacobian(q_bar, qs))'*D*inv(jacobian(q_bar, qs));
    
    % modified variables
    qb1 = sym('qb1', 'real');
    qb2 = sym('qb2', 'real');
    qN = sym('qN', 'real');
    q_bar = [qb1; qb2; qN];
    q = H\(q_bar - b);
    V_bar = subs(V_tot, [q1 q2 q3], q');
    D_bar = subs(D_bar, [q1 q2 q3], q');
end