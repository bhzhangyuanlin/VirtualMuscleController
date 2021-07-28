clear all;
       p11 =      -743.3;
       p21 =        2170;
       p31 =       -2472;
       p41 =        1362;
       p51 =        -357;
       p61 =       36.87;
       p1 =  -2.886e+04;
       p2 =   8.491e+04;
       p3 =  -9.678e+04;
       p4 =   5.321e+04;
       p5 =  -1.406e+04;
       p6 =   1442;
counter=0;
f=zeros(1,51);
a=zeros(1,51);
t=zeros(1,51);
LineW=1.3;
for x=0.3:0.01:0.8
      counter=counter+1;
   f(counter)= p11*x^5 + p21*x^4 + p31*x^3 + p41*x^2 + p51*x + p61;
  a(counter) = p1*x^5 + p2*x^4 + p3*x^3 + p4*x^2 + p5*x + p6;
  t(counter)=x;
  end;
  figure(10)


[ax,h1,h2]=plotyy(t,f,t,a,'plot');
set(get(ax(1),'Ylabel'),'string','f(Hz)','color','b') %y1
set(get(ax(2),'Ylabel'),'string','A(mm)','color','r') %y2
xlabel('V(m/s)');
set(h1,'linestyle','-','color','b','LineWidth',LineW);   
set(h2,'linestyle','- -','color','r','LineWidth',LineW);
legend([h1 h2],'f','A') %标注两条线
set(ax(1),'Ycolor','b') %设定Y轴的颜色为蓝色
set(ax(2),'Ycolor','r') %设定Y轴的颜色为红色
set(ax(1),'ytick',[1:0.2:3]); %设置y轴间隔
set(ax(2),'ytick',[10:2:30])
set(gca,'xtick',[0.3:0.05:0.8]);% 设置x轴范围
grid on;

