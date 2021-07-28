
function outyita=roboticfish401(f,Kp_1,Kp_2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               说明                                       %

%   在原来实验的基础上，发现了运动速度和频率以及鱼身体刚度的关系（用三次函数来 % 
%拟合）。修改程序在原来fish_normal.m文件基础上增加了速度变化，以及用频率的线性调
%节来改变机器鱼的速度以达到对指定速度的跟随。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               变量声明                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%   运动学部分变量声明   %%%%%%%%%%%%%%%%%%%%%%%%%%%

PHI          = pi/2;%1.34;         %第一关节第二关节相位差


theta_10f    = 0.3144;       %第一关节相对于鱼体头部的转动角度幅值
theta_21f    = 1.0123;       %第二关节相对于第一关节转动的角度幅值

V_c1c        = 0;            %鱼体身体部分质心的垂直于鱼尾表面的线速度分量
V_c1         = 0;            %鱼体身体部分质心处线速度
V_c1_v       = [0,0];        %鱼体身体部分质心处线速度，矢量
V_1c        = 0;            %鱼体身体部分质心与水流的合速度垂直于鱼尾表面的线速度分量
V_1         = 0;            %鱼体身体部分质心处线速度与水流的合速度
V_1_v        = [0,0];       %鱼体身体部分质心处线速度与水流的合速度，矢量



V_c2c        = 0;            %鱼体尾部质心的垂直于鱼尾表面的线速度分量
V_c2         = 0;            %鱼体尾部质心处线速度
V_c2_v       = [0,0];        %鱼体尾部质心处线速度，矢量
V_2c        = 0;            %鱼体尾部质心与水流的合速度垂直于鱼尾表面的线速度分量
V_2         = 0;            %鱼体尾部质心处线速度与水流的合速度
V_2_v       = [0,0];        %鱼体尾部质心处线速度与水流的合速度，矢量



theta_10     = 0;            %鱼体身体相对于鱼体头部角度
theta_21     = theta_21f * sin(PHI);
                             %鱼体尾部相对于鱼体身体角度
theta_20     = theta_21f * sin(PHI);
                             %鱼体尾部相对于鱼体头部的角度

theta_10_1o_g= 0;            %给定的theta_1_0对t的一阶导数
theta_21_1o_g= 0;            %给定的theta_2_1对t的一阶导数
theta_20_1o_g= 0;            %给定的theta_2_0对t的一阶导数

X_0          = 0;
X_0_1o       = 0;            %鱼体头部横坐标的一阶导数即行进速度
U            = X_0_1o;       %鱼体行进速度
X_1          = 0.32 ;%鱼身体中心横坐标X_1         
X_2          = 0.52 ;%鱼尾中心横坐标X_2
pos_1        = [0.32,0];%鱼身体中心位置，矢量
pos_2        = [0.52,0];%鱼尾中心位置，矢量theta_10_1o  = 0;            %theta_1_0对t的一阶导数

theta_10_1o=0;
theta_21_1o  = 0;            %theta_2_1对t的一阶导数
theta_20_1o  = 0;            %theta_2_0对t的一阶导数


theta_10_g   = 0;            %给定的theta_1_0对t的一阶导数
theta_21_g   = 0;            %给定的theta_2_1对t的一阶导数
theta_20_g   = 0;            %给定的theta_2_0对t的一阶导数

theta_10_2o  = 0;            %theta_1_0对t的二阶导数
theta_21_2o  = 0;            %theta_2_1对t的二阶导数
theta_20_2o  = 0;            %theta_2_0对t的二阶导数
X_0_2o       = 0;            %鱼体头部横坐标的二阶导数即行进加速度

I_2_v        = [0;0];        %鱼体尾部单位法向量


%%%%%%%%%%%%   动力学部分变量声明   %%%%%%%%%%%%%%

rou          = 1000;         %流体密度
C_0          = 0.9;          %鱼体头部阻力系数
C_1          = 1.54;         %鱼体身体部分阻力系数
C_2          = 1.54;         %鱼体尾部阻力系数

Q_1          = 0;            %相对于广义坐标theta_10的非保守力
Q_2          = 0;            %相对于广义坐标theta_21的非保守力
Q_3          = 0;            %相对于广义坐标X_0的非保守力

F_0          = 0;            %鱼体头部水动力
F_0_v        = [0,0];        %鱼体头部水动力，矢量 
F_1_s        = 0;            %鱼体身体部分水动力，标量
F_1_v        = [0,0];        %鱼体身体部分水动力，矢量 
F_1x         = 0;            %鱼体身体部分水动力在X轴上的分量
F_2_s        = 0;            %鱼体尾部水动力,标量
F_2_v        = [0,0];        %鱼体尾部水动力，矢量 
F_2x         = 0;            %鱼体尾部水动力在X轴上的分量

M_1          = 0;            %第一个关节的转矩
M_2          = 0;            %第二个关节的转矩

P_1          = 0;             %第一个关节的的功率
P_2          = 0;             %第二个关节的功率
P_usefu      = 0;

W_total      = 0;            %鱼体游动时消耗的总功
W_useful     = 0;            %鱼体游动时的有用功
yita         = 0;            %机器鱼的推进效率

%%%%%%%%%%   拉格朗日方程系数变量声明   %%%%%%%%%%%

D_11         = 0;            %拉格朗日系数D_11
D_12         = 0;            %拉格朗日系数D_12
D_13         = 0;            %拉格朗日系数D_13
D_21         = 0;            %拉格朗日系数D_21
D_22         = 0;            %拉格朗日系数D_22
D_23         = 0;            %拉格朗日系数D_23
D_31         = 0;            %拉格朗日系数D_31
D_32         = 0;            %拉格朗日系数D_32
D_33         = 0;            %拉格朗日系数D_33

D_111        = 0;            %拉格朗日系数D_111
D_122        = 0;            %拉格朗日系数D_122
D_211        = 0;            %拉格朗日系数D_211
D_222        = 0;            %拉格朗日系数D_222
D_311        = 0;            %拉格朗日系数D_311
D_322        = 0;            %拉格朗日系数D_322

D_112        = 0;            %拉格朗日系数D_112
D_212        = 0;            %拉格朗日系数D_212
D_312        = 0;            %拉格朗日系数D_312

DD_1         = zeros(3);     %3 X 3矩阵
DD_1_inv     = zeros(3);     %DD_1矩阵的逆
DD_2         = zeros(3,2);   %3 X 2矩阵
DD_3         = zeros(3,1);   %3 X 1矩阵

%%%%%%%%%%%%   鱼体身体相关变量声明   %%%%%%%%%%%%%

m_0          = 4.5;          %鱼体头部质量
m_1          = 0.4;          %鱼体身体部分质量
m_2          = 1;            %鱼体尾部质量

A_0x         = 0.03;         %鱼体头部与速度方向垂直面积
A_1          = 0.05;         %身体部分投影面积
A_2          = 0.05;         %尾部投影面积

L_0          = 0.2;          %鱼体头部长度

L_1          = 0.24;         %鱼体身体部分长度
L_c1         = 0.12;         %鱼体身体部分质心距离第一个关节的长度
L_d1         = 0.12;         %鱼体身体部分水动力作用点距离第一个关节的长度

L_2          = 0.16;         %鱼体尾部长度
L_c2         = 0.08;         %鱼体尾部质心距离第二个关节的长度
L_d2         = 0.08;         %鱼体尾部水动力作用点距离第二个关节的长度


J_c1         = m_1 * L_1 * L_1 / 12;
                             %鱼体身体部分以质心为参考点的转动惯量
J_c2         = m_2 * L_2 * L_2 / 12;
                             %鱼体尾部以质心为参考点的转动惯量

%%%%%%%%%%%%     其他变量声明     %%%%%%%%%%%%%%

time         = 30;           %仿真时间
dt           = 0.01;         %仿真时间间隔
times        = 1 / dt;       %时间与此系数相乘即得到数组个数

counter      = 0;            %计算循环次数，将数据存入到相应数组中的元素
x_axis       = zeros(1,time * times);
                             %画图用X轴坐标，与时间轴上取点个数相同
                             
                             
                             
F_1_s_t      = zeros(1,time * times);
                             %身体部分受力大小，标量，无方向                             
F_1_x_t      = zeros(1,time * times);
F_1_y_t      = zeros(1,time * times);
V_c1c_t      = zeros(1,time * times);
                             %身体部分摆动产生的前进推力与其中心的线速度，与时间轴上取点个数相同                             
V_1c_t      = zeros(1,time * times);                             
                              %鱼体中心的线速度与水流速度的合速度，与时间轴上取点个数相同
                             
                             
F_2_s_t      = zeros(1,time * times);
                             %尾部受力大小，标量，无方向
F_2x_t       = zeros(1,time * times);
F_2y_t       = zeros(1,time * times);
V_c2c_t      = zeros(1,time * times);
                             %尾部摆动产生的前进推力与其中心的线速度，与时间轴上取点个数相同
V_2c_t      = zeros(1,time * times);                             
                              %鱼尾中心的线速度与水流速度的合速度，与时间轴上取点个数相同                             
                             

theta_10_t   = zeros(1,time * times);
theta_20_t   = zeros(1,time * times);
theta_21_t   = zeros(1,time * times);
X_0_t        = zeros(1,time * times);
                             %记录各个关节的角度与鱼的位移，与时间轴上取点个数相同

theta_10_g_t = zeros(1,time * times);
theta_21_g_t = zeros(1,time * times);
theta_20_g_t = zeros(1,time * times);
                             %记录给定信号的取值与时间轴上取点个数相同

theta_10_1o_t= zeros(1,time * times);
theta_21_1o_t= zeros(1,time * times);
theta_20_1o_t= zeros(1,time * times);
X_0_1o_t     = zeros(1,time * times);
                             %记录各个关节的角速度与鱼前进速度，与时间轴上取点个数相同

theta_10_2o_t= zeros(1,time * times);
theta_21_2o_t= zeros(1,time * times);
theta_20_2o_t= zeros(1,time * times);
X_0_2o_t     = zeros(1,time * times);
                             %记录各个关节的角加速度与鱼的前进加速度，与时间轴上取点个数相同

M_1_t        = zeros(1,time * times);
M_2_t        = zeros(1,time * times);
                             %记录两个关节的转矩，与时间轴上取点个数相同
                             
P_1_t          = zeros(1,time * times);
P_2_t          = zeros(1,time * times);
                             %记录两个关节的功率，与时间轴上取点个数相同
                             
P_useful_t          = zeros(1,time * times);


temp         = zeros(3,1);   %用于计算拉格朗日方程的解的时候的临时变量，因为MATLAB中赋值操作左侧的数组中不允许出现“；”

start_time   = 0;           %用于计算效率的开始与结束时间,计算两个周期的总功与有用功
end_time     = start_time + 30;

%%%%%%%%%%%%%   PID参数变量声明   %%%%%%%%%%%%%%

%Kp_1         = 5;          %第一个关节（身体）的PID参数
Ki_1         = 0;
Kd_1         = 0;
err_1        = 0;
d_err_1      = 0;
sumerr_1     = 0;
last_err_1   = 0;

%Kp_2         = 10;          %第二个关节（尾部）的PID参数
Ki_2         = 0;
Kd_2         = 0;
err_2        = 0;
d_err_2      = 0;
sumerr_2     = 0;
last_err_2   = 0;


v_g=0;                              %给定速度 
f_g=0;                              %给定频率（由给定速度算出）
Kp_g=0;                             %给定刚度（由给定速度算出）
v_g_t=zeros(1,time*times);          %给定速度数组
f_g_t=zeros(1,time*times);          %给定频率数组
f_t=zeros(1,time*times);
k=1;                                %给定速度的计数器
       
counter_a=0;                        %求平均速度的计数开始节点
counter_b=0;                        %求平均速度的计数结束节点
j=1;                                %循环变量
sum=0;                              %一段时间内的速度求和
v_av=0;                             %一段时间内的平均速度
t_v_b=0;
t_v_a=0;
v_x(1)=0;
v_y(1)=0;
%%%%%求平均实时功率的变量%%%%%%%%%%%%%
P_a=0;
P_b=0;                                                                  %%%%计数器末端
W_tot=0;                                                               %%%%某一段时间内的总功
W_use=0;                                                              %%%%某一段时间内的有用功
yita_k=0;                                                                %%%%暂存的效率变量
yita_t=zeros(1,time*times);                                     %%%%实时效率
i=1;                                                                       %%%循环变量
t_b=0;
t_a=0;
yita_x(1)=0;
yita_y(1)=0;

%%闭环频率%%
k_u=0.005;
Kp_g_t=zeros(1,time*times);
Kp_1_t=zeros(1,time*times);
Kp_2_t=zeros(1,time*times);
Kp_1_c=0.1;
Kp_2_c=0.1;
t_p=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               运算开始                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = dt : dt : time                                                     %从0.01秒到10秒
    counter =counter + 1;                                                  %计数器加1          
    if(t==dt)                                                              %if语句实现的功能为在开始时给定鱼一个速度并算出相应的频率和第一关节的刚度（速度、频率、刚度关系已知）
          v_g=0.4;
         % f_g=round(21.48*v_g^3-21.64*v_g^2+13.37*v_g-0.6503);
         % f=1.8;  
      %  %  Kp_g=200;
   end; 
       if(floor(t)==15*k)                                                   %这里if语句实现功能为每过5s给定速度变化一次。floor函数为向下取整函数，解决t速度精度累积误差问题。
            v_g=v_g+0.1;
            t_p=t;
            k=k+1;
         end;   

    Kp_g_t(counter)=Kp_g;
    f_t(counter)=f;
    v_g_t(counter)=v_g;
    f_g_t(counter)=f_g;                            
    x_axis(counter)=t;   
    k_g_t(counter)=k;
    
    if(k==1)
      if(t-t_p>3)
        Kp_1=Kp_1+(v_g-U)*Kp_1_c;
        Kp_2=Kp_2+(v_g-U)*Kp_2_c;
     end;

    else
         Kp_1=Kp_1+(v_g-U)*Kp_1_c;
        Kp_2=Kp_2+(v_g-U)*Kp_2_c;
end;
  Kp_1_t(counter)=Kp_1;
  Kp_2_t(counter)=Kp_2;
    
    theta_10_g = theta_10f * sin(2 * pi * f * t);                          %theta_10_1o随时间变化函数
    theta_21_g = theta_21f * sin(2 * pi * f * t + PHI);                    %theta_21_1o随时间变化的函数
    theta_20_g = theta_10_g + theta_21_g;                                  %theta_20_1o随时间变化的函数
    
    theta_10_g_t(counter) = theta_10_g;
    theta_21_g_t(counter) = theta_21_g;
    theta_20_g_t(counter) = theta_20_g;
    
    err_1 = theta_10_g - theta_10;                                         %第一个关节当前角速度与设定角速度之间差值
    d_err_1 = err_1 - last_err_1;                                          %第一个关节当前偏差与上次偏差的差值
    last_err_1 = err_1;                                                    %第一个关节将当前角速度偏差设定为上次偏差

  
    M_1 = Kp_1 * err_1 + Kd_1 * d_err_1;                                   %第一个关节的转矩
    M_1_t(counter) = M_1;
    
    err_2 = theta_21_g - theta_21;                                         %第二个关节当前角速度与设定角速度之间差值
    d_err_2 = err_2 - last_err_2;                                          %第二个关节当前偏差与上次偏差的差值
    last_err_2 = err_2;                                                    %第二个关节将当前角速度偏差设定为上次偏差
    
    M_2 = Kp_2 * err_2 + Kd_2 * d_err_2;                                   %第二个关节的转矩
    M_2_t(counter) = M_2; 
    
    
    P_1 = M_1 * theta_10_1o ;                                              %第一个关节的功率  
    P_2 = M_2 * theta_21_1o;                                               %第二个关节的功率
    
    P_1_t(counter) = P_1;
    P_2_t(counter) = P_2;
    
    D_11 = J_c1 + J_c2 + m_1 * L_c1 * L_c1 + m_2 * (L_1 * L_1 + L_c2 * L_c2 + 2 * L_1 * L_c2 * cos(theta_21));
    D_12 = J_c2 + m_2 * (L_c2 * L_c2 + L_1 * L_c2 * cos (theta_21));
    D_21 = D_12;
    D_13 = -m_1 * L_c1 * sin(theta_10) - m_2 * (L_1 * sin(theta_10) + L_c2 * sin(theta_20));
    D_31 = D_13;
    D_22 = m_2 * L_c2 * L_c2 + J_c2;
    D_23 = -m_2 * L_c2 * sin(theta_20);
    D_32 = D_23;
    D_33 = m_0 + m_1 + m_2;
    D_111 = 0;
    D_222 = 0;
    D_122 = -m_2 * L_1 * L_c2 * sin(theta_21);
    D_211 = D_122;
    D_311 = -m_1 * L_c1 * cos(theta_10) - m_2 * (L_1 * cos(theta_10) + L_c2 * cos(theta_20));
    D_322 = -m_2 * L_c2 * cos(theta_20);
    D_112 = -2 * m_2 * L_1 * L_c2 * sin(theta_21);
    D_212 = 0;
    D_312 = -2 * m_2 * L_c2 * cos(theta_20);
    
    DD_1 = [D_11 D_12 D_13;
            D_21 D_22 D_23;
            D_31 D_32 D_33];
%     DD_1_inv = inv(DD_1);
    
    DD_2 = [D_111 D_122;
            D_211 D_222;
            D_311 D_322];
        
    DD_3 = [D_112;
            D_212;
            D_312];
        
             
    A = csvread('wave_position.csv');
   m=size(A,1);
    for i=1:m
       B(i,1)=A(i,1);
       B(i,2)=A(i,2);
    end

    for i=1:m
       D1(i)=(B(i,1)-pos_1(1)).^2+(B(i,2)-pos_1(2)).^2;
       D2(i)=(B(i,1)-pos_2(1)).^2+(B(i,2)-pos_2(2)).^2;
    end
   [d1,m1]=min(D1);
   [d2,m2]=min(D2);
   X1(1)=A(m1,3);
   X1(2)=A(m1,4);
   X2(1)=A(m2,3);
   X2(2)=A(m2,4);
      
             
  V_c1_v = [X_0_1o - L_1 * sin(theta_10) * theta_10_1o ,L_1 * cos(theta_10) * theta_10_1o ];
                                                                           %V_c1的横纵坐标
    I_1 = [sin(theta_10);-cos(theta_10)];                                  %鱼体身体部分速度的单位法向量
    V_c1c = V_c1_v * I_1;                                                  %鱼体身体部分法向速度
    V_c1c_t(counter) = V_c1c;
    
    
    
    V_1_v = V_c1_v + X1;
    V_1c = V_1_v * I_1;
    V_1c_t(counter) = V_1c;
    
    
    
    
    
    
    
    
    F_1_s = 0.5 * 0.5 * rou * C_1 * V_1c * V_1c;                         %鱼体身体部分受力，标量，为正
    F_1_s_t(counter) = F_1_s;
    if(V_1c > 0)
        F_1_v = [-0.5 * rou * C_1 * V_1c * V_1c * A_1 * sin(theta_10),0.5 * rou * C_1 * V_1c * V_1c * A_1 * cos(theta_10)];
                                                                           %当摆动方向和法向量方向相同时，身体部分产生的水动力
        F_1_x_t(counter) = F_1_v(1);
        F_1x = F_1_v(1);
        F_1_y_t(counter) = F_1_v(2);
    else
        F_1_v = [0.5 * rou * C_1 * V_1c * V_1c * A_1 * sin(theta_10),-0.5 * rou * C_1 * V_1c * V_1c * A_1 * cos(theta_10)];
                                                                           %当摆动方向和法向量方向相反时，身体部分产生的水动力
        F_1_x_t(counter) = F_1_v(1);
        F_1x = F_1_v(1);
        F_1_y_t(counter) = F_1_v(2);
    end
    
    
    V_c2_v = [X_0_1o - L_1 * sin(theta_10) * theta_10_1o - L_c2 * sin(theta_20) * theta_20_1o,L_1 * cos(theta_10) * theta_10_1o + L_c2 * cos(theta_20) * theta_20_1o];
                                                                           %V_c2的横纵坐标
    I_2 = [sin(theta_20);-cos(theta_20)];                                  %鱼体尾部速度的单位法向量
    V_c2c = V_c2_v * I_2;                                                  %鱼体尾部法向速度
    V_c2c_t(counter) = V_c2c;
    
    
    
        
    V_2_v = V_c2_v + X2;
    V_2c = V_2_v * I_2;
    V_2c_t(counter) = V_2c;
    
    
    F_2_s = 0.5 * 0.5 * rou * C_2 * V_2c * V_2c;                         %鱼体尾部受力，标量，为正
    F_2_s_t(counter) = F_2_s;
    if(V_2c > 0)
        F_2_v = [-0.5 * rou * C_2 * V_2c * V_2c * A_2 * sin(theta_20),0.5 * rou * C_2 * V_2c * V_2c * A_2 * cos(theta_20)];
                                                                           %当摆动方向和法向量方向相同时，尾部产生的水动力
        F_2x_t(counter) = F_2_v(1);
        F_2x = F_2_v(1);
        F_2y_t(counter) = F_2_v(2);
    else
        F_2_v = [0.5 * rou * C_2 * V_2c * V_2c * A_2 * sin(theta_20),-0.5 * rou * C_2 * V_2c * V_2c * A_2 * cos(theta_20)];
                                                                           %当摆动方向和法向量方向相反时，尾部产生的水动力
        F_2x_t(counter) = F_2_v(1);
        F_2x = F_2_v(1);
        F_2y_t(counter) = F_2_v(2);

    end
    F_0 = -0.5 * rou * C_0 * X_0_1o * X_0_1o * A_0x;
    F_0_v = [F_0,0];
    
    Q_1 = F_2_v * [-L_1 * sin(theta_10) - L_d2 * sin(theta_20);L_1 * cos(theta_10) + L_d2 * cos(theta_20)] + M_1;
    Q_2 = F_2_v * [-L_d2 * sin(theta_20);L_d2 * cos(theta_20)] + M_2;
    Q_3 = F_0_v * [1;0] + F_2_v * [1;0];
        
    temp = DD_1 \ ([Q_1;Q_2;Q_3] - DD_2 * [theta_10_1o ^ 2;theta_21_1o ^ 2] - DD_3 * (theta_10_1o * theta_21_1o));      
    
    pos_1         = pos_1 + V_c1_v * dt;%鱼身体中心位置，矢量
     pos_2         = pos_2 + V_c2_v * dt;%鱼身体中心位置，矢量


    theta_10_2o = temp(1,1);
    theta_21_2o = temp(2,1);
    X_0_2o = temp(3,1);
    
    theta_10_2o_t(counter) = theta_10_2o;
    theta_21_2o_t(counter) = theta_21_2o;
    X_0_2o_t(counter) = X_0_2o;
    
    theta_10_1o = theta_10_1o + theta_10_2o * dt;
    theta_21_1o = theta_21_1o + theta_21_2o * dt;
    X_0_1o = X_0_1o + X_0_2o * dt; 
    U = X_0_1o;
    
    theta_10_1o_t(counter) = theta_10_1o;
    theta_21_1o_t(counter) = theta_21_1o;
    theta_20_1o_t(counter) = theta_21_1o + theta_10_1o;
    X_0_1o_t(counter) = X_0_1o;
  
    theta_10 = theta_10 + theta_10_1o * dt;
    theta_10_t(counter) = theta_10;
    
    theta_21 = theta_21 + theta_21_1o * dt;
    theta_21_t(counter) = theta_21; 

    
    theta_20 = theta_10 + theta_21;
    theta_20_t(counter) = theta_20;
    X_0 = X_0 + X_0_1o * dt; 
    X_0_t(counter) = X_0;
    P_useful = ( F_1x + F_2x ) * U;
    P_useful_t(counter) = P_useful;
    
    if t >= start_time && t <= end_time
           if(P_1>0) && (P_2>0)
        W_total = W_total + (P_1+P_2) * dt;
        else if (P_1>0) && (P_2<0)
              W_total = W_total+ (P_1) * dt;  
            else if (P_1<0) && (P_2>0)
                    W_total = W_total + (P_2) * dt;  
                else W_total = W_total ;
                end;
            end;
           end;
        W_useful = W_useful + P_useful * dt;
    end
    
    %%%%平均功率实时处理，计算每一间隔时间的平均效率，并做出效率时间图
    %%%和速度的变化放在一块看
    
    t_b=t;
       if(P_1>0) && (P_2>0)
        W_tot = W_tot + (P_1+P_2) * dt;
        else if (P_1>0) && (P_2<0)
              W_tot = W_tot + (P_1) * dt;  
            else if (P_1<0) && (P_2>0)
                    W_tot = W_tot + (P_2) * dt;  
                else W_tot = W_tot ;
                end;
            end;
       end;

    W_use=W_use+P_useful*dt;
    if(round(4/f,2)-t_b+t_a<=1e-10)                                   %%算平均功率的时间为50*0.01=0.5s，大约为一个周期，也可以是其他的数。
       yita_k=W_use/W_tot;
       i=i+1;
       yita_y(i)=yita_k;
     
       yita_x(i)=t_b;
      t_b_t(i)=t_b;
      t_a_t(i)=t_a;
      W_tot=0;
      W_use=0;
      t_a=t_b;
     end;  
                                                                           %求平均速度的目的：处理数据，使数据曲线看起来更平缓，对机器鱼性能没有影响，可以省去。
    sum=sum+X_0_1o;                                                        %速度累加
    t_v_b=t;                                                     %求平均速度的计数结束节点
    if(round(2/f,2)-t_v_b+t_v_a<=1e-10)                                                      %每计数100次，即每过（100*0.01=1s)求一次平均速度
       counter_b=t_v_b*1/dt;
       counter_a=t_v_a*1/dt;
        v_av=sum/(counter_b-counter_a);
        j=j+1;
        v_y(j)=v_av;
        v_x(j)=t_v_b;
        sum=0;
        t_v_a=t_v_b;
    end;
    %counter=counter_a~counter_b这段内的速度相同（重新赋值为平均速度）
                                                                                                                                                             
    end;
    
            


yita = W_useful / W_total;
disp(['f = ',num2str(f)]);
disp(['Kp_1 = ',num2str(Kp_1)]);
disp(['Kp_2 = ',num2str(Kp_2)]);
disp(['W_total = ',num2str(W_total)]);
disp(['W_useful = ',num2str(W_useful)]);
disp(['yita = ',num2str(yita)]);

figure (1);
plot(x_axis,v_g_t,'r',v_x,v_y,'b--'); %鱼体头部横坐标的一阶导数即行进速度
xlabel('t(time/s)');
ylabel('speed(m/s)');
legend('V_{gt}','U_{t}');
axis([0 30 0 0.7]);
grid on;

figure(2);
plot(x_axis,P_1_t,'r',x_axis,P_2_t,'b--');
xlabel('t(time/s)');
ylabel('P(w)');
legend('P_{1t}','P_{2t}');
axis([0 15 -15 15]);
grid on;

figure(3);
plot(x_axis,theta_10_g_t,'r',x_axis,theta_10_t,'b--');
legend('\theta_{1g}','\theta_{1}');
xlabel('t(time/s)');
ylabel('\theta(rad)');

%axis([0 15 -0.25 0.25]);
grid on;

figure(4);
plot(x_axis,theta_21_g_t,'r-',x_axis,theta_21_t,'b--');
legend('\theta_{2g}','\theta_{2}');
xlabel('t(time/s)');
ylabel('\theta(rad)');
axis([0 15 -1.5 1.5]);
grid on;

figure(5);
plot(x_axis,theta_10_1o_t,'r-',x_axis,theta_21_1o_t,'b--');
xlabel('t(time/s)');
ylabel('theta_ dot');
legend('theta1_ dot','theta2_ dot');
grid on;
%axis([0 15 -10 10]);

figure(6);
plot(x_axis,M_1_t,'r-',x_axis,M_2_t,'b--');
xlabel('time/s');
ylabel('M_t');
legend('M_{1t}','M_{2t}');
grid on;

figure(7);
plot(yita_x,yita_y);
xlabel('t(time/s)');
ylabel('\eta');
grid on;
outyita=yita;




% 
% figure (1);                                                                  
% plot(x_axis,V_c1c_t); %鱼体身体部分质心的垂直于鱼尾表面的线速度分量
% xlabel('time/s')
% ylabel('m/s')
% axis([0 30 -0.3 0.3]);
% title('The line speed  vertical component of fish body centroid');
% 
% figure (2);                                                                  
% plot(x_axis,F_1_x_t);%鱼体身体部分水动力在X轴上的分量
% xlabel('time/s')
% ylabel('N')
% axis([0 time -0.15 0.15]);
% title('hydrodynamic component of fish body  in the X-axis');
% 
% figure (3);
% plot(x_axis,V_c2c_t); %鱼体尾部质心的垂直于鱼尾表面的线速度分量
% xlabel('time/s')
% ylabel('m/s')
% axis([0 30 -0.3 0.3]);
% title('The line speed  vertical component of fishtail centroid');

% figure (4);
% plot(x_axis,F_2x_t);%鱼体尾部水动力在X轴上的分量
% xlabel('time/s')
% ylabel('N')
% axis([25 30 -0.5 0.5]);
% legend('F2x')
% title('hydrodynamic component of fishtail  in the X-axis');
% % 
% figure (5);
% plot(x_axis,X_0_t); %鱼体头部横坐标
% xlabel('time/s')
% ylabel('m')
% axis([0 time -3 3]);
% title('abscissa of fish head');
% % 

% % 
% figure (7);
% plot(x_axis,theta_20_t);  %鱼体尾部相对于鱼体头部的角度
% xlabel('time/s')
% ylabel('rad')
% axis([0 time -0.8 0.8]);
% title('angle of fishtail relative to fish head');
% 
% figure (8);
% plot(x_axis,theta_21_g_t,x_axis,theta_21_t,':','LineWidth',2);%给定的鱼体身体相对于鱼体头部角度
% legend('θ','θ1')
% xlabel('time/s')
% ylabel('rad')
% axis([10 15 -0.8 0.8]);

% plot();  %鱼体身体相对于鱼体头部角度
% xlabel('time/s')
% ylabel('rad')
% axis([0 time -0.2 0.2]);
% title('angle of fish body relative to fish head');
% 
% figure (9);
% plot(x_axis,theta_21_t);%鱼体尾部相对于鱼体身体角度
% xlabel('time/s')
% ylabel('rad')
% axis([0 time -0.8 0.8]);
% title('angle of fishtail relative to fish body');
% 
% figure (10);

% plot(x_axis,theta_10_1o_t);%鱼体身体相对于鱼体头部角速度
% xlabel('time/s')
% ylabel('rad/s')
% axis([0 time -1 1]);
% title('angular velocity of fish body relative to fish head');
% 
% figure (11);
% plot(x_axis,theta_21_1o_t);%鱼体尾部相对于鱼体身体角速度
% xlabel('time/s')
% ylabel('rad/s')
% axis([0 time -7 7]);
% title('angular velocity of fishtail relative to fish body');
% 
% figure (12)
% plot(x_axis,theta_10_2o_t); %鱼体尾部相对于鱼体头部角加速度
% xlabel('time/s')
% ylabel('rad/s^2')
% axis([0 time -10 10]);
% title('angular acceleration of fishtail relative to fish head');
% 
% figure (13)
% plot(x_axis,theta_21_2o_t);%鱼体尾部相对于鱼体身体角加速度
% xlabel('time/s')
% ylabel('rad/s^2')
% axis([0 time -25 25]);
% title('angular acceleration of fishtail relative to fish head');
% 
% figure (14)
% plot(x_axis,X_0_2o_t); %鱼体行进加速度
% xlabel('time/s')
% ylabel('m/s^2')
% axis([0 time -0.2 0.2]);
% title('traveling acceleration of fish');
% 
% figure (15)
% plot(x_axis,theta_10_g_t);%给定的鱼体身体相对于鱼体头部角度
% xlabel('time/s')
% ylabel('rad')
% axis([0 time -0.2 0.2]);
% title('given angle of fish body relative to fish head');
% 
% figure (16)
% plot(x_axis,theta_21_g_t);%给定的鱼体尾部相对于鱼体身体角度
% xlabel('time/s')
% ylabel('rad')
% axis([0 time -1 1]);
% title('given angle of fishtail relative to fish body');
% 
% figure (17)
% plot(x_axis,theta_20_g_t);%给定的鱼体尾部相对于鱼体头部角度
% xlabel('time/s')
% ylabel('rad')
% axis([0 time -1 1]);
% title('given angle of fishtail relative to fish body');
% 
% figure (18)
% plot(x_axis,M_1_t);%第一个关节的转矩
% xlabel('time/s')
% ylabel('N*m')
% axis([0 time -0.4 0.4]);
% title('torque of the first joint');
% 
% figure (19)
% plot(x_axis,M_2_t);%第二个关节的转矩
% xlabel('time/s')
% ylabel('N*m')
% axis([0 time -0.2 0.2]);
% title('torque of the second join');
% 
% figure (20)
% plot(x_axis,F_1_s_t);%鱼体身体部分受力，标量，为正
% xlabel('time/s')
% ylabel('N')
% axis([0 time -10 25]);
% title('fish body force, scalar, positive');
% 
% figure (21)
% plot(x_axis,F_2_s_t);%鱼体尾部受力，标量，为正
% xlabel('time/s')
% ylabel('N')
% axis([0 time -10 25]);
% title('fishtail force, scalar, positive');
% 
% figure (22)
% plot(x_axis,M_1_t,x_axis,theta_10_1o_t);%%,'r',x_axis,x_zeros); %第一关节转矩及鱼身体部分相对于鱼头部的角速度
% legend('M_1','\theta_{10}dot')
% xlabel('time/s')
% ylabel('M/N*m w/rad^-1')
% axis([0 time -1.2 1.2]);
% title('torque of the first joint and angular velocity of fish body relative to fish head');
% 
% 
% figure (23)
% plot(x_axis,M_2_t,x_axis,theta_20_1o_t);%%,'r',x_axis,x_zeros); %第二关节转矩及鱼尾部相对于鱼身体的角速度
% legend('M_2','\theta_20dot')
% xlabel('time/s')
% ylabel('M/N*m w/rad*s^-1')
% axis([0 time -4 4]);
% title('torque of the second joint and angular velocity of fishtail relative to fish boay');
% 
% 
% 
% figure (24)
% plot(x_axis,M_1_t,x_axis,theta_10_t);%%,'r',x_axis,x_zeros); %第一关节转矩及鱼身体部分相对于鱼头部的角度
% legend('M_1','\theta_{10}')
% xlabel('time/s')
% ylabel('M/N*m θ/rad*s^-1')
% axis([0 time -0.3 0.3]);

% 
% 
% 
% figure (25)
% plot(x_axis,M_2_t,x_axis,theta_21_t);%%,'r',x_axis,x_zeros);  %第二关节转矩及鱼尾部相对于鱼身体部分的角度
% legend('M_2','\theta_{21}')
% xlabel('time/s')
% ylabel('M/N*m θ/rad*s^-1')
% axis([0 time -0.7 0.7]);
% % title('torque of the first joint and angule of fishtail relative to fish body');



