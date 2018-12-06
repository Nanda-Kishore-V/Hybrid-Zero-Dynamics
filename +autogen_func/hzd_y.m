function y = hzd_y(q1,q2,q3,dq1,dq2,dq3)
%HZD_Y
%    Y = HZD_Y(Q1,Q2,Q3,DQ1,DQ2,DQ3)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    29-Nov-2018 18:19:32

t2 = q1.^2;
t3 = t2.^2;
y = [q1.*4.09989700010116e-2-q3+t2.*3.708933218415968e-1+t3.*1.313290086870686e-1+pi-q1.*t2.*2.112513582720696+q1.*t3.*1.018662501650811e1-t2.*t3.*2.093315400356285e1+q1.*t2.*t3.*1.141188754230027e1-2.628405923403671;q1.*1.504461068584838+q2-q3+t2.*3.564855763392445-t3.*2.30674622682309+pi-q1.*t2.*6.94550413273528+q1.*t3.*2.092553787789513e1-t2.*t3.*4.431576949045961e1+q1.*t2.*t3.*1.836055100659802e1-2.977019839728937];
