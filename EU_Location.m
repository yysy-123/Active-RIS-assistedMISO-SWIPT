function value=EU_Location(x0,y0,R,K)
% 输入依次是：圆心横纵坐标，半径和点的数量
theta=0:0.001:360;
% 利用极坐标得到圆的坐标
Circle1=x0+R*cos(theta);
Circle2=y0+R*sin(theta);
% 画圆
% plot(Circle1,Circle2,'r')
% 随机生成num_Dian个半径
r=R*sqrt(rand(1,K));
% 得到生成点的角度，并利用极坐标形式画出点
seta=2*pi*rand(1,K);
% 得到点的坐标
x=x0+r.*cos(seta);
y=y0+r.*sin(seta);
value=[x.',y.'];