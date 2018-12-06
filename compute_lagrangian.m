function [] = compute_lagrangian()
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
    
    save('+data/potential_energy.mat', 'V_tot');
    save('+data/kinetic_energy.mat', 'K_tot')
end