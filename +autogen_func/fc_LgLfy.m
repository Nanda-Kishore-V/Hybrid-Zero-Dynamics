function LgLfy = fc_LgLfy(q1,q2,q3,dq1,dq2,dq3)
%FC_LGLFY
%    LGLFY = FC_LGLFY(Q1,Q2,Q3,DQ1,DQ2,DQ3)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    27-Nov-2018 16:53:27

t2 = q1.*2.0;
t3 = q2.*-2.0+t2;
t4 = cos(t3);
t5 = q1-q3;
t6 = cos(t5);
t7 = q1.^2;
t9 = q2.*-1.0;
t10 = q1+t9;
t8 = cos(t10);
t11 = t8.^2;
t12 = q3.*-1.0;
t13 = q1+t12;
t14 = cos(t13);
t15 = t11.*4.0;
t16 = t14.^2;
t17 = t16.*8.0;
t18 = t15+t17-2.5e1;
t19 = 1.0./t18;
t20 = t7.^2;
LgLfy = reshape([-(q1.*7.0e1-t4.*1.0e3+t6.*2.146e3-t7.*2.457e3+q1.*t6.*1.4e2-t6.*t7.*4.914e3+1.1573e4)./(t4.*2.5e3+cos(q3.*-2.0+t2).*5.0e3-2.375e4),t19.*(t14.*2.0+1.0).*(q1.*3.519489393747758e4-t7.*5.699593430458233e4+t8.*1.28e4-t20.*6.048e4-q1.*t7.*7.9616e4+9.617491034755131e3).*1.25e-4,t19.*(t11.*-4.0+t8.*t14.*8.0+2.5e1).*-4.0e-1-t19.*(t8+t14).*(q1.*(7.0./1.0e2)-t7.*2.457+7.3e1./1.0e3).*1.6,(t16.*-3.2e1+t8.*t14.*1.6e1+1.0e2)./(t11.*2.0e1+t16.*4.0e1-1.25e2)-t19.*(t6+cos(q1-q2)).*(q1.*-8.798723484369395+t7.*1.424898357614558e1+t20.*1.512e1+q1.*t7.*1.9904e1-2.404372758688783)],[2,2]);
