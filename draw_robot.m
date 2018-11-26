function [] = draw_robot(qe)
% constants
l = 0.5;
r = 1.0;
M_T = 10.0;
M_H = 15.0;
m = 5.0;
g = 9.81;

% state
q1 = -qe(1);
q2 = -qe(2);
q3 = -qe(3);
p0h = qe(4);
p0v = qe(5);

% unit vectors
i = [1; 0];
j = [0; 1];
e1 = cos(q1)*j - sin(q1)*i;
e2 = -cos(q2)*j + sin(q2)*i;
e3 = cos(q3)*j - sin(q3)*i;

% link ends
P0 = p0h*i + p0v*j;
P1 = P0 + r*e1;
P2 = P1 + r*e2;
P3 = P1 + l*e3;

% link center of masses
G1 = P0 + r*e1/2;
G2 = P1 + r*e2/2;

% bounds
max_height = l+r;
xbnds = max_height*[-1.2, 5.2];
ybnds = max_height*[-0.2, 1.2];

% Colors:
color_ground = [118,62,12]/255;
color_stance = [200,60,60]/255;
color_swing = [60,60,200]/255;
color_torso = [160, 80, 160]/255;

hold off;

plot(xbnds, [0,0],'LineWidth',6,'Color',color_ground);

hold on;

% plot links
plot([P0(1), P1(1)], [P0(2), P1(2)],'LineWidth',4,'Color',color_stance);
plot([P1(1), P2(1)], [P1(2), P2(2)],'LineWidth',4,'Color',color_swing);
plot([P1(1), P3(1)], [P1(2), P3(2)],'LineWidth',4,'Color',color_torso);

% plot joints
plot(P0(1), P0(2),'k.','MarkerSize',30);
plot(P1(1), P1(2),'k.','MarkerSize',30);
plot(P2(1), P2(2),'k.','MarkerSize',30);
plot(P3(1), P3(2),'k.','MarkerSize',30);

% plot COM
plot(G1(1), G1(2),'ko','MarkerSize',8,'LineWidth',2);
plot(G2(1), G2(2),'ko','MarkerSize',8,'LineWidth',2);

axis([xbnds,ybnds]);
axis equal;
axis off;

end