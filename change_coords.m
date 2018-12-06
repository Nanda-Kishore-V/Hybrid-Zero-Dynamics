function [D_body, V_body] = change_coords(H, b)
    % variables
    q1 = sym('q1', 'real');
    q2 = sym('q2', 'real');
    q3 = sym('q3', 'real');

    dq1 = sym('dq1', 'real');
    dq2 = sym('dq2', 'real');
    dq3 = sym('dq3', 'real');
    
    qs = [q1; q2; q3];
    
    load('+data/D.mat', 'D');
    load('+data/potential_energy.mat', 'V_tot');
    q_bar = H*qs + b;
    D_bar = inv(jacobian(q_bar, qs))'*D*inv(jacobian(q_bar, qs));
    
    % modified variables
    qb1 = sym('qb1', 'real');
    qb2 = sym('qb2', 'real');
    qN = sym('qN', 'real');
    q_bar = [qb1; qb2; qN];
    q = H\(q_bar - b);
    V_body = subs(V_tot, [q1 q2 q3], q');
    D_body = subs(D_bar, [q1 q2 q3], q');
    
    save('+data/V_body.mat', 'V_body');
    save('+data/D_body.mat', 'D_body');
end