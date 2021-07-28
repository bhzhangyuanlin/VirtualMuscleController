       clear;
       p1 =   4.918e-05  ;
       p2 =   -0.002163 ;
       p3 =     0.01528  ;
       p4 =      0.3096 ;
       counter=0;
       kp1x=zeros(1,30);
       t=zeros(1,30);
       for x=1:1:30
           counter=counter+1;
           kp1x(counter)= p1*x^3 + p2*x^2 + p3*x + p4;
           t(counter)=x;
       end;
figure(11)
plot(t,kp1x,'k'); %鱼体头部横坐标的一阶导数即行进速度
xlabel('K_1');
ylabel('\eta');
legend('K_2=16');
grid on;
axis([ 0 30 0.1 0.35]);

       a1 =      0.4546  ;
       b1 =       0.106 ;
       c1 =     -0.4912 ;
       a2 =      0.1437 ;
       b2 =      0.2576 ;
       c2 =     -0.2129  ;
       a3 =     0.02943  ;
       b3 =      0.4385 ; 
       c3 =     -0.5094  ;
       a4 =    0.007552 ;
       b4 =      0.7245  ;
       c4 =       3.706  ;
       counter=0;
       kp2x=zeros(1,30);
       t=zeros(1,59);
         for x=1:0.5:30
           counter=counter+1;
           kp2x(counter)= a1*sin(b1*x+c1) + a2*sin(b2*x+c2) + a3*sin(b3*x+c3) + a4*sin(b4*x+c4);
           t(counter)=x;
       end;
       
figure(12)
plot(t,kp2x,'k'); %鱼体头部横坐标的一阶导数即行进速度
xlabel('K_2');
ylabel('\eta');
legend('K_1=5');
grid on;
axis([ 8 30 0.3 0.36]);


       a11 =      0.4819  ;
       b11 =      0.1423  ;
       c11=       1.475  ;
       a21 =      0.1649  ;
       b21 =      0.2645 ;
       c21 =      -1.791 ;
       a31 =     0.02543;  
       b31 =      0.6651 ; 
       c31 =      0.8482 ;
       a41 =      0.0152  ;
       b41 =      0.9401  ;
       c41 =      -2.811  ;
       a51 =    0.006833 ;
       b51 =       1.116  ;
       c51 =     -0.1565  ;
       counter=0;
     H1x=zeros(1,41);
       t=zeros(1,41);
        for x=-10:0.5:10
           counter=counter+1;
           H1x(counter)= a11*sin(b11*x+c11) + a21*sin(b21*x+c21) + a31*sin(b31*x+c31) + a41*sin(b41*x+c41)+a51*sin(b51*x+c51);
           t(counter)=x;
       end;
       
figure(13)
plot(t,H1x,'k'); %鱼体头部横坐标的一阶导数即行进速度
xlabel('H_1');
ylabel('\eta');
legend('H_2=9');
grid on;

    

       p11 =  -9.922e-07 ;
       p21 =   2.773e-06  ;
       p31 =   0.0001412 ;
       p41 =  -0.0007167 ;
       p51 =    0.001589  ;
       p61 =       0.316 ;
       
     counter=0;
     H2x=zeros(1,41);
       t=zeros(1,41);
        for x=-10:0.5:10
           counter=counter+1;
           H2x(counter)= p11*x^5 + p21*x^4 + p31*x^3 + p41*x^2 + p51*x + p61;
           t(counter)=x;
       end;
figure(14)
plot(t,H2x,'k'); %鱼体头部横坐标的一阶导数即行进速度
xlabel('H_2');
ylabel('\eta');
legend('H_1=1');
axis([-10 10 0.2 0.36]);
grid on;
       
       
           
