function value=EU_Location(x0,y0,R,K)
% ���������ǣ�Բ�ĺ������꣬�뾶�͵������
theta=0:0.001:360;
% ���ü�����õ�Բ������
Circle1=x0+R*cos(theta);
Circle2=y0+R*sin(theta);
% ��Բ
% plot(Circle1,Circle2,'r')
% �������num_Dian���뾶
r=R*sqrt(rand(1,K));
% �õ����ɵ�ĽǶȣ������ü�������ʽ������
seta=2*pi*rand(1,K);
% �õ��������
x=x0+r.*cos(seta);
y=y0+r.*sin(seta);
value=[x.',y.'];