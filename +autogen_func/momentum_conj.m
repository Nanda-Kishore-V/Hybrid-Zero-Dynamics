function out1 = momentum_conj(qb1,qb2,qN,dqb1,dqb2,dqN)
%MOMENTUM_CONJ
%    OUT1 = MOMENTUM_CONJ(QB1,QB2,QN,DQB1,DQB2,DQN)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    27-Nov-2018 18:09:00

t2 = qb2.*-1.0;
t3 = qb1+t2;
t4 = cos(t3);
t5 = cos(qb1);
t6 = t4.*1.25;
t7 = t4.*2.5;
t8 = t5.*5.0;
out1 = dqb1.*(t5.*2.5+t6-1.875)-dqN.*(t7+t8-1.75e1).*1.0+dqb1.*(t7+t8-3.75).*5.0e-1-dqN.*(t4.*5.0+t5.*1.0e1-3.5e1).*5.0e-1-dqb2.*(t6-6.25e-1).*1.0-dqb2.*(t7-1.25).*5.0e-1;