function Lf2y = fc_Lf2y(q1,q2,q3,dq1,dq2,dq3)
%FC_LF2Y
%    LF2Y = FC_LF2Y(Q1,Q2,Q3,DQ1,DQ2,DQ3)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    27-Nov-2018 16:53:28

t2 = q1.*2.0;
t9 = q2.*-2.0;
t5 = -t9;
t3 = -t5;
t42 = q1.*2.0;
t25 = -t42;
t14 = -t25;
t11 = -t14;
t4 = -t11;
t6 = q3.*-1.0;
t7 = q3.*-2.0;
t8 = dq1.^2;
t10 = dq2.^2;
t12 = t2+t7;
t13 = cos(t12);
t15 = dq3.^2;
t16 = q1-q3;
t17 = sin(t16);
t18 = t9+t14;
t19 = t7+t14;
t20 = sin(t18);
t21 = q1-q2;
t22 = sin(t21);
t23 = q2-q3;
t24 = sin(t23);
t26 = pi.^2;
t27 = q1.^2;
t28 = t2+t9;
t29 = cos(t28);
t30 = t29.*2.0;
t31 = t13.*4.0;
t32 = t30+t31-1.9e1;
t33 = 1.0./t32;
t34 = q1+t9;
t35 = sin(t34);
t36 = t35.*9.81e2;
t37 = q1+t7;
t38 = sin(t37);
t39 = t38.*1.962e3;
t40 = sin(q1);
t41 = t40.*9.81e3;
t43 = t10.*t22.*1.0e2;
t44 = t15.*t17.*-2.0e2;
Lf2y = [t8.*(q1.*4.914-7.0./1.0e2)-((sin(q3+t3+t4).*-9.81e2-sin(t4+t6).*1.1772e4+sin(t5+t6).*9.81e2+sin(q3).*1.0791e4+t8.*t17.*2.3e3+t10.*t24.*1.0e2+t15.*sin(t4+t7).*2.0e2+t8.*sin(q1+q3+t9).*2.0e2+t10.*sin(q2+q3+t11).*1.0e2).*1.0)./(t13.*2.0e2+cos(t2+t3).*1.0e2-9.5e2)+t33.*(q1.*1.4e-3-t27.*4.914e-2+1.46e-3).*(t36+t39+t41+t43+t44-t8.*t20.*1.0e2-t8.*sin(t19).*2.0e2);t8.*(q1.*6.2592e4-t26.*3.11e2+t27.*1.19424e5-q1.*t26.*5.67e2+q1.*t27.*1.2096e5-1.4528e4).*(-3.125e-4)-((sin(q2+t7+t14).*1.962e3+sin(q2.*-1.0+t14).*1.0791e4+sin(q2+t7).*1.962e3-sin(q2).*9.81e3-t8.*t22.*2.1e3+t10.*t20.*1.0e2-t15.*t24.*2.0e2-t8.*sin(q1+q2+t7).*4.0e2+t15.*sin(q2+q3+t25).*2.0e2).*1.0)./(cos(t18).*1.0e2+cos(t19).*2.0e2-9.5e2)+t33.*(q1.*-1.099840435546174e-1+t27.*1.781122947018198e-1+q1.*t27.*2.488e-1+t27.^2.*1.89e-1-3.005465948360978e-2).*(t36+t39+t41+t43+t44-t8.*sin(t7+t42).*2.0e2-t8.*sin(t9+t42).*1.0e2)];
