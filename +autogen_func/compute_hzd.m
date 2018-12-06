function [] = compute_hzd()
% unfinished
    precomputation()

    import autogen_func.swing_model
    import autogen_func.impact_model
    import autogen_func.hzd_LgLfy
    import autogen_func.hzd_Lf2y

    H = [1 0 -1; 0 1 -1; 1 0 0];
    b = [pi; pi; 0];

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

    matlabFunction(jacobian(K, dqN), 'File', '+autogen_func/momentum_conj.m', 'vars', [qb1, qb2, qN, dqb1, dqb2, dqN]);

    theta = qN;
    hd1 = poly2sym(hd(1, :), qN);
    hd2 = poly2sym(hd(2, :), qN);
    h1 = qb1 - hd1;
    h2 = qb2 - hd2;
    h = [h1; h2];
    k1 = jacobian(theta, qb)*inv([jacobian(h, qb); gamma_0])*[0; 0; 1];
    matlabFunction(k1, 'File', '+autogen_func/k1', 'vars', [qb1, qb2, qN]);
    k1 = subs(k1, [qb1, qb2], [hd1, hd2]);
    k2 = -jacobian(V_body, qN);
    matlabFunction(k2, 'File', '+autogen_func/k2', 'vars', [qb1, qb2, qN]);
end