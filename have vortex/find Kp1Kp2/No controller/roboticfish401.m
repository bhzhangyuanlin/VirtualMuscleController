
function outyita=roboticfish401(f,Kp_1,Kp_2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               ˵��                                       %

%   ��ԭ��ʵ��Ļ����ϣ��������˶��ٶȺ�Ƶ���Լ�������նȵĹ�ϵ�������κ����� % 
%��ϣ����޸ĳ�����ԭ��fish_normal.m�ļ��������������ٶȱ仯���Լ���Ƶ�ʵ����Ե�
%�����ı��������ٶ��Դﵽ��ָ���ٶȵĸ��档
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               ��������                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%   �˶�ѧ���ֱ�������   %%%%%%%%%%%%%%%%%%%%%%%%%%%

PHI          = pi/2;%1.34;         %��һ�ؽڵڶ��ؽ���λ��


theta_10f    = 0.3144;       %��һ�ؽ����������ͷ����ת���Ƕȷ�ֵ
theta_21f    = 1.0123;       %�ڶ��ؽ�����ڵ�һ�ؽ�ת���ĽǶȷ�ֵ

V_c1c        = 0;            %�������岿�����ĵĴ�ֱ����β��������ٶȷ���
V_c1         = 0;            %�������岿�����Ĵ����ٶ�
V_c1_v       = [0,0];        %�������岿�����Ĵ����ٶȣ�ʸ��
V_1c        = 0;            %�������岿��������ˮ���ĺ��ٶȴ�ֱ����β��������ٶȷ���
V_1         = 0;            %�������岿�����Ĵ����ٶ���ˮ���ĺ��ٶ�
V_1_v        = [0,0];       %�������岿�����Ĵ����ٶ���ˮ���ĺ��ٶȣ�ʸ��



V_c2c        = 0;            %����β�����ĵĴ�ֱ����β��������ٶȷ���
V_c2         = 0;            %����β�����Ĵ����ٶ�
V_c2_v       = [0,0];        %����β�����Ĵ����ٶȣ�ʸ��
V_2c        = 0;            %����β��������ˮ���ĺ��ٶȴ�ֱ����β��������ٶȷ���
V_2         = 0;            %����β�����Ĵ����ٶ���ˮ���ĺ��ٶ�
V_2_v       = [0,0];        %����β�����Ĵ����ٶ���ˮ���ĺ��ٶȣ�ʸ��



theta_10     = 0;            %�����������������ͷ���Ƕ�
theta_21     = theta_21f * sin(PHI);
                             %����β���������������Ƕ�
theta_20     = theta_21f * sin(PHI);
                             %����β�����������ͷ���ĽǶ�

theta_10_1o_g= 0;            %������theta_1_0��t��һ�׵���
theta_21_1o_g= 0;            %������theta_2_1��t��һ�׵���
theta_20_1o_g= 0;            %������theta_2_0��t��һ�׵���

X_0          = 0;
X_0_1o       = 0;            %����ͷ���������һ�׵������н��ٶ�
U            = X_0_1o;       %�����н��ٶ�
X_1          = 0.32 ;%���������ĺ�����X_1         
X_2          = 0.52 ;%��β���ĺ�����X_2
pos_1        = [0.32,0];%����������λ�ã�ʸ��
pos_2        = [0.52,0];%��β����λ�ã�ʸ��theta_10_1o  = 0;            %theta_1_0��t��һ�׵���

theta_10_1o=0;
theta_21_1o  = 0;            %theta_2_1��t��һ�׵���
theta_20_1o  = 0;            %theta_2_0��t��һ�׵���


theta_10_g   = 0;            %������theta_1_0��t��һ�׵���
theta_21_g   = 0;            %������theta_2_1��t��һ�׵���
theta_20_g   = 0;            %������theta_2_0��t��һ�׵���

theta_10_2o  = 0;            %theta_1_0��t�Ķ��׵���
theta_21_2o  = 0;            %theta_2_1��t�Ķ��׵���
theta_20_2o  = 0;            %theta_2_0��t�Ķ��׵���
X_0_2o       = 0;            %����ͷ��������Ķ��׵������н����ٶ�

I_2_v        = [0;0];        %����β����λ������


%%%%%%%%%%%%   ����ѧ���ֱ�������   %%%%%%%%%%%%%%

rou          = 1000;         %�����ܶ�
C_0          = 0.9;          %����ͷ������ϵ��
C_1          = 1.54;         %�������岿������ϵ��
C_2          = 1.54;         %����β������ϵ��

Q_1          = 0;            %����ڹ�������theta_10�ķǱ�����
Q_2          = 0;            %����ڹ�������theta_21�ķǱ�����
Q_3          = 0;            %����ڹ�������X_0�ķǱ�����

F_0          = 0;            %����ͷ��ˮ����
F_0_v        = [0,0];        %����ͷ��ˮ������ʸ�� 
F_1_s        = 0;            %�������岿��ˮ����������
F_1_v        = [0,0];        %�������岿��ˮ������ʸ�� 
F_1x         = 0;            %�������岿��ˮ������X���ϵķ���
F_2_s        = 0;            %����β��ˮ����,����
F_2_v        = [0,0];        %����β��ˮ������ʸ�� 
F_2x         = 0;            %����β��ˮ������X���ϵķ���

M_1          = 0;            %��һ���ؽڵ�ת��
M_2          = 0;            %�ڶ����ؽڵ�ת��

P_1          = 0;             %��һ���ؽڵĵĹ���
P_2          = 0;             %�ڶ����ؽڵĹ���
P_usefu      = 0;

W_total      = 0;            %�����ζ�ʱ���ĵ��ܹ�
W_useful     = 0;            %�����ζ�ʱ�����ù�
yita         = 0;            %��������ƽ�Ч��

%%%%%%%%%%   �������շ���ϵ����������   %%%%%%%%%%%

D_11         = 0;            %��������ϵ��D_11
D_12         = 0;            %��������ϵ��D_12
D_13         = 0;            %��������ϵ��D_13
D_21         = 0;            %��������ϵ��D_21
D_22         = 0;            %��������ϵ��D_22
D_23         = 0;            %��������ϵ��D_23
D_31         = 0;            %��������ϵ��D_31
D_32         = 0;            %��������ϵ��D_32
D_33         = 0;            %��������ϵ��D_33

D_111        = 0;            %��������ϵ��D_111
D_122        = 0;            %��������ϵ��D_122
D_211        = 0;            %��������ϵ��D_211
D_222        = 0;            %��������ϵ��D_222
D_311        = 0;            %��������ϵ��D_311
D_322        = 0;            %��������ϵ��D_322

D_112        = 0;            %��������ϵ��D_112
D_212        = 0;            %��������ϵ��D_212
D_312        = 0;            %��������ϵ��D_312

DD_1         = zeros(3);     %3 X 3����
DD_1_inv     = zeros(3);     %DD_1�������
DD_2         = zeros(3,2);   %3 X 2����
DD_3         = zeros(3,1);   %3 X 1����

%%%%%%%%%%%%   ����������ر�������   %%%%%%%%%%%%%

m_0          = 4.5;          %����ͷ������
m_1          = 0.4;          %�������岿������
m_2          = 1;            %����β������

A_0x         = 0.03;         %����ͷ�����ٶȷ���ֱ���
A_1          = 0.05;         %���岿��ͶӰ���
A_2          = 0.05;         %β��ͶӰ���

L_0          = 0.2;          %����ͷ������

L_1          = 0.24;         %�������岿�ֳ���
L_c1         = 0.12;         %�������岿�����ľ����һ���ؽڵĳ���
L_d1         = 0.12;         %�������岿��ˮ�������õ�����һ���ؽڵĳ���

L_2          = 0.16;         %����β������
L_c2         = 0.08;         %����β�����ľ���ڶ����ؽڵĳ���
L_d2         = 0.08;         %����β��ˮ�������õ����ڶ����ؽڵĳ���


J_c1         = m_1 * L_1 * L_1 / 12;
                             %�������岿��������Ϊ�ο����ת������
J_c2         = m_2 * L_2 * L_2 / 12;
                             %����β��������Ϊ�ο����ת������

%%%%%%%%%%%%     ������������     %%%%%%%%%%%%%%

time         = 30;           %����ʱ��
dt           = 0.01;         %����ʱ����
times        = 1 / dt;       %ʱ�����ϵ����˼��õ��������

counter      = 0;            %����ѭ�������������ݴ��뵽��Ӧ�����е�Ԫ��
x_axis       = zeros(1,time * times);
                             %��ͼ��X�����꣬��ʱ������ȡ�������ͬ
                             
                             
                             
F_1_s_t      = zeros(1,time * times);
                             %���岿��������С���������޷���                             
F_1_x_t      = zeros(1,time * times);
F_1_y_t      = zeros(1,time * times);
V_c1c_t      = zeros(1,time * times);
                             %���岿�ְڶ�������ǰ�������������ĵ����ٶȣ���ʱ������ȡ�������ͬ                             
V_1c_t      = zeros(1,time * times);                             
                              %�������ĵ����ٶ���ˮ���ٶȵĺ��ٶȣ���ʱ������ȡ�������ͬ
                             
                             
F_2_s_t      = zeros(1,time * times);
                             %β��������С���������޷���
F_2x_t       = zeros(1,time * times);
F_2y_t       = zeros(1,time * times);
V_c2c_t      = zeros(1,time * times);
                             %β���ڶ�������ǰ�������������ĵ����ٶȣ���ʱ������ȡ�������ͬ
V_2c_t      = zeros(1,time * times);                             
                              %��β���ĵ����ٶ���ˮ���ٶȵĺ��ٶȣ���ʱ������ȡ�������ͬ                             
                             

theta_10_t   = zeros(1,time * times);
theta_20_t   = zeros(1,time * times);
theta_21_t   = zeros(1,time * times);
X_0_t        = zeros(1,time * times);
                             %��¼�����ؽڵĽǶ������λ�ƣ���ʱ������ȡ�������ͬ

theta_10_g_t = zeros(1,time * times);
theta_21_g_t = zeros(1,time * times);
theta_20_g_t = zeros(1,time * times);
                             %��¼�����źŵ�ȡֵ��ʱ������ȡ�������ͬ

theta_10_1o_t= zeros(1,time * times);
theta_21_1o_t= zeros(1,time * times);
theta_20_1o_t= zeros(1,time * times);
X_0_1o_t     = zeros(1,time * times);
                             %��¼�����ؽڵĽ��ٶ�����ǰ���ٶȣ���ʱ������ȡ�������ͬ

theta_10_2o_t= zeros(1,time * times);
theta_21_2o_t= zeros(1,time * times);
theta_20_2o_t= zeros(1,time * times);
X_0_2o_t     = zeros(1,time * times);
                             %��¼�����ؽڵĽǼ��ٶ������ǰ�����ٶȣ���ʱ������ȡ�������ͬ

M_1_t        = zeros(1,time * times);
M_2_t        = zeros(1,time * times);
                             %��¼�����ؽڵ�ת�أ���ʱ������ȡ�������ͬ
                             
P_1_t          = zeros(1,time * times);
P_2_t          = zeros(1,time * times);
                             %��¼�����ؽڵĹ��ʣ���ʱ������ȡ�������ͬ
                             
P_useful_t          = zeros(1,time * times);


temp         = zeros(3,1);   %���ڼ����������շ��̵Ľ��ʱ�����ʱ��������ΪMATLAB�и�ֵ�������������в�������֡�����

start_time   = 0;           %���ڼ���Ч�ʵĿ�ʼ�����ʱ��,�����������ڵ��ܹ������ù�
end_time     = start_time + 30;

%%%%%%%%%%%%%   PID������������   %%%%%%%%%%%%%%

%Kp_1         = 5;          %��һ���ؽڣ����壩��PID����
Ki_1         = 0;
Kd_1         = 0;
err_1        = 0;
d_err_1      = 0;
sumerr_1     = 0;
last_err_1   = 0;

%Kp_2         = 10;          %�ڶ����ؽڣ�β������PID����
Ki_2         = 0;
Kd_2         = 0;
err_2        = 0;
d_err_2      = 0;
sumerr_2     = 0;
last_err_2   = 0;


v_g=0;                              %�����ٶ� 
f_g=0;                              %����Ƶ�ʣ��ɸ����ٶ������
Kp_g=0;                             %�����նȣ��ɸ����ٶ������
v_g_t=zeros(1,time*times);          %�����ٶ�����
f_g_t=zeros(1,time*times);          %����Ƶ������
f_t=zeros(1,time*times);
k=1;                                %�����ٶȵļ�����
       
counter_a=0;                        %��ƽ���ٶȵļ�����ʼ�ڵ�
counter_b=0;                        %��ƽ���ٶȵļ��������ڵ�
j=1;                                %ѭ������
sum=0;                              %һ��ʱ���ڵ��ٶ����
v_av=0;                             %һ��ʱ���ڵ�ƽ���ٶ�
t_v_b=0;
t_v_a=0;
v_x(1)=0;
v_y(1)=0;
%%%%%��ƽ��ʵʱ���ʵı���%%%%%%%%%%%%%
P_a=0;
P_b=0;                                                                  %%%%������ĩ��
W_tot=0;                                                               %%%%ĳһ��ʱ���ڵ��ܹ�
W_use=0;                                                              %%%%ĳһ��ʱ���ڵ����ù�
yita_k=0;                                                                %%%%�ݴ��Ч�ʱ���
yita_t=zeros(1,time*times);                                     %%%%ʵʱЧ��
i=1;                                                                       %%%ѭ������
t_b=0;
t_a=0;
yita_x(1)=0;
yita_y(1)=0;

%%�ջ�Ƶ��%%
k_u=0.005;
Kp_g_t=zeros(1,time*times);
Kp_1_t=zeros(1,time*times);
Kp_2_t=zeros(1,time*times);
Kp_1_c=0.1;
Kp_2_c=0.1;
t_p=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               ���㿪ʼ                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = dt : dt : time                                                     %��0.01�뵽10��
    counter =counter + 1;                                                  %��������1          
    if(t==dt)                                                              %if���ʵ�ֵĹ���Ϊ�ڿ�ʼʱ������һ���ٶȲ������Ӧ��Ƶ�ʺ͵�һ�ؽڵĸնȣ��ٶȡ�Ƶ�ʡ��նȹ�ϵ��֪��
          v_g=0.4;
         % f_g=round(21.48*v_g^3-21.64*v_g^2+13.37*v_g-0.6503);
         % f=1.8;  
      %  %  Kp_g=200;
   end; 
       if(floor(t)==15*k)                                                   %����if���ʵ�ֹ���Ϊÿ��5s�����ٶȱ仯һ�Ρ�floor����Ϊ����ȡ�����������t�ٶȾ����ۻ�������⡣
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
    
    theta_10_g = theta_10f * sin(2 * pi * f * t);                          %theta_10_1o��ʱ��仯����
    theta_21_g = theta_21f * sin(2 * pi * f * t + PHI);                    %theta_21_1o��ʱ��仯�ĺ���
    theta_20_g = theta_10_g + theta_21_g;                                  %theta_20_1o��ʱ��仯�ĺ���
    
    theta_10_g_t(counter) = theta_10_g;
    theta_21_g_t(counter) = theta_21_g;
    theta_20_g_t(counter) = theta_20_g;
    
    err_1 = theta_10_g - theta_10;                                         %��һ���ؽڵ�ǰ���ٶ����趨���ٶ�֮���ֵ
    d_err_1 = err_1 - last_err_1;                                          %��һ���ؽڵ�ǰƫ�����ϴ�ƫ��Ĳ�ֵ
    last_err_1 = err_1;                                                    %��һ���ؽڽ���ǰ���ٶ�ƫ���趨Ϊ�ϴ�ƫ��

  
    M_1 = Kp_1 * err_1 + Kd_1 * d_err_1;                                   %��һ���ؽڵ�ת��
    M_1_t(counter) = M_1;
    
    err_2 = theta_21_g - theta_21;                                         %�ڶ����ؽڵ�ǰ���ٶ����趨���ٶ�֮���ֵ
    d_err_2 = err_2 - last_err_2;                                          %�ڶ����ؽڵ�ǰƫ�����ϴ�ƫ��Ĳ�ֵ
    last_err_2 = err_2;                                                    %�ڶ����ؽڽ���ǰ���ٶ�ƫ���趨Ϊ�ϴ�ƫ��
    
    M_2 = Kp_2 * err_2 + Kd_2 * d_err_2;                                   %�ڶ����ؽڵ�ת��
    M_2_t(counter) = M_2; 
    
    
    P_1 = M_1 * theta_10_1o ;                                              %��һ���ؽڵĹ���  
    P_2 = M_2 * theta_21_1o;                                               %�ڶ����ؽڵĹ���
    
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
                                                                           %V_c1�ĺ�������
    I_1 = [sin(theta_10);-cos(theta_10)];                                  %�������岿���ٶȵĵ�λ������
    V_c1c = V_c1_v * I_1;                                                  %�������岿�ַ����ٶ�
    V_c1c_t(counter) = V_c1c;
    
    
    
    V_1_v = V_c1_v + X1;
    V_1c = V_1_v * I_1;
    V_1c_t(counter) = V_1c;
    
    
    
    
    
    
    
    
    F_1_s = 0.5 * 0.5 * rou * C_1 * V_1c * V_1c;                         %�������岿��������������Ϊ��
    F_1_s_t(counter) = F_1_s;
    if(V_1c > 0)
        F_1_v = [-0.5 * rou * C_1 * V_1c * V_1c * A_1 * sin(theta_10),0.5 * rou * C_1 * V_1c * V_1c * A_1 * cos(theta_10)];
                                                                           %���ڶ�����ͷ�����������ͬʱ�����岿�ֲ�����ˮ����
        F_1_x_t(counter) = F_1_v(1);
        F_1x = F_1_v(1);
        F_1_y_t(counter) = F_1_v(2);
    else
        F_1_v = [0.5 * rou * C_1 * V_1c * V_1c * A_1 * sin(theta_10),-0.5 * rou * C_1 * V_1c * V_1c * A_1 * cos(theta_10)];
                                                                           %���ڶ�����ͷ����������෴ʱ�����岿�ֲ�����ˮ����
        F_1_x_t(counter) = F_1_v(1);
        F_1x = F_1_v(1);
        F_1_y_t(counter) = F_1_v(2);
    end
    
    
    V_c2_v = [X_0_1o - L_1 * sin(theta_10) * theta_10_1o - L_c2 * sin(theta_20) * theta_20_1o,L_1 * cos(theta_10) * theta_10_1o + L_c2 * cos(theta_20) * theta_20_1o];
                                                                           %V_c2�ĺ�������
    I_2 = [sin(theta_20);-cos(theta_20)];                                  %����β���ٶȵĵ�λ������
    V_c2c = V_c2_v * I_2;                                                  %����β�������ٶ�
    V_c2c_t(counter) = V_c2c;
    
    
    
        
    V_2_v = V_c2_v + X2;
    V_2c = V_2_v * I_2;
    V_2c_t(counter) = V_2c;
    
    
    F_2_s = 0.5 * 0.5 * rou * C_2 * V_2c * V_2c;                         %����β��������������Ϊ��
    F_2_s_t(counter) = F_2_s;
    if(V_2c > 0)
        F_2_v = [-0.5 * rou * C_2 * V_2c * V_2c * A_2 * sin(theta_20),0.5 * rou * C_2 * V_2c * V_2c * A_2 * cos(theta_20)];
                                                                           %���ڶ�����ͷ�����������ͬʱ��β��������ˮ����
        F_2x_t(counter) = F_2_v(1);
        F_2x = F_2_v(1);
        F_2y_t(counter) = F_2_v(2);
    else
        F_2_v = [0.5 * rou * C_2 * V_2c * V_2c * A_2 * sin(theta_20),-0.5 * rou * C_2 * V_2c * V_2c * A_2 * cos(theta_20)];
                                                                           %���ڶ�����ͷ����������෴ʱ��β��������ˮ����
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
    
    pos_1         = pos_1 + V_c1_v * dt;%����������λ�ã�ʸ��
     pos_2         = pos_2 + V_c2_v * dt;%����������λ�ã�ʸ��


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
    
    %%%%ƽ������ʵʱ��������ÿһ���ʱ���ƽ��Ч�ʣ�������Ч��ʱ��ͼ
    %%%���ٶȵı仯����һ�鿴
    
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
    if(round(4/f,2)-t_b+t_a<=1e-10)                                   %%��ƽ�����ʵ�ʱ��Ϊ50*0.01=0.5s����ԼΪһ�����ڣ�Ҳ����������������
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
                                                                           %��ƽ���ٶȵ�Ŀ�ģ��������ݣ�ʹ�������߿�������ƽ�����Ի���������û��Ӱ�죬����ʡȥ��
    sum=sum+X_0_1o;                                                        %�ٶ��ۼ�
    t_v_b=t;                                                     %��ƽ���ٶȵļ��������ڵ�
    if(round(2/f,2)-t_v_b+t_v_a<=1e-10)                                                      %ÿ����100�Σ���ÿ����100*0.01=1s)��һ��ƽ���ٶ�
       counter_b=t_v_b*1/dt;
       counter_a=t_v_a*1/dt;
        v_av=sum/(counter_b-counter_a);
        j=j+1;
        v_y(j)=v_av;
        v_x(j)=t_v_b;
        sum=0;
        t_v_a=t_v_b;
    end;
    %counter=counter_a~counter_b����ڵ��ٶ���ͬ�����¸�ֵΪƽ���ٶȣ�
                                                                                                                                                             
    end;
    
            


yita = W_useful / W_total;
disp(['f = ',num2str(f)]);
disp(['Kp_1 = ',num2str(Kp_1)]);
disp(['Kp_2 = ',num2str(Kp_2)]);
disp(['W_total = ',num2str(W_total)]);
disp(['W_useful = ',num2str(W_useful)]);
disp(['yita = ',num2str(yita)]);

figure (1);
plot(x_axis,v_g_t,'r',v_x,v_y,'b--'); %����ͷ���������һ�׵������н��ٶ�
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
% plot(x_axis,V_c1c_t); %�������岿�����ĵĴ�ֱ����β��������ٶȷ���
% xlabel('time/s')
% ylabel('m/s')
% axis([0 30 -0.3 0.3]);
% title('The line speed  vertical component of fish body centroid');
% 
% figure (2);                                                                  
% plot(x_axis,F_1_x_t);%�������岿��ˮ������X���ϵķ���
% xlabel('time/s')
% ylabel('N')
% axis([0 time -0.15 0.15]);
% title('hydrodynamic component of fish body  in the X-axis');
% 
% figure (3);
% plot(x_axis,V_c2c_t); %����β�����ĵĴ�ֱ����β��������ٶȷ���
% xlabel('time/s')
% ylabel('m/s')
% axis([0 30 -0.3 0.3]);
% title('The line speed  vertical component of fishtail centroid');

% figure (4);
% plot(x_axis,F_2x_t);%����β��ˮ������X���ϵķ���
% xlabel('time/s')
% ylabel('N')
% axis([25 30 -0.5 0.5]);
% legend('F2x')
% title('hydrodynamic component of fishtail  in the X-axis');
% % 
% figure (5);
% plot(x_axis,X_0_t); %����ͷ��������
% xlabel('time/s')
% ylabel('m')
% axis([0 time -3 3]);
% title('abscissa of fish head');
% % 

% % 
% figure (7);
% plot(x_axis,theta_20_t);  %����β�����������ͷ���ĽǶ�
% xlabel('time/s')
% ylabel('rad')
% axis([0 time -0.8 0.8]);
% title('angle of fishtail relative to fish head');
% 
% figure (8);
% plot(x_axis,theta_21_g_t,x_axis,theta_21_t,':','LineWidth',2);%�����������������������ͷ���Ƕ�
% legend('��','��1')
% xlabel('time/s')
% ylabel('rad')
% axis([10 15 -0.8 0.8]);

% plot();  %�����������������ͷ���Ƕ�
% xlabel('time/s')
% ylabel('rad')
% axis([0 time -0.2 0.2]);
% title('angle of fish body relative to fish head');
% 
% figure (9);
% plot(x_axis,theta_21_t);%����β���������������Ƕ�
% xlabel('time/s')
% ylabel('rad')
% axis([0 time -0.8 0.8]);
% title('angle of fishtail relative to fish body');
% 
% figure (10);

% plot(x_axis,theta_10_1o_t);%�����������������ͷ�����ٶ�
% xlabel('time/s')
% ylabel('rad/s')
% axis([0 time -1 1]);
% title('angular velocity of fish body relative to fish head');
% 
% figure (11);
% plot(x_axis,theta_21_1o_t);%����β�����������������ٶ�
% xlabel('time/s')
% ylabel('rad/s')
% axis([0 time -7 7]);
% title('angular velocity of fishtail relative to fish body');
% 
% figure (12)
% plot(x_axis,theta_10_2o_t); %����β�����������ͷ���Ǽ��ٶ�
% xlabel('time/s')
% ylabel('rad/s^2')
% axis([0 time -10 10]);
% title('angular acceleration of fishtail relative to fish head');
% 
% figure (13)
% plot(x_axis,theta_21_2o_t);%����β���������������Ǽ��ٶ�
% xlabel('time/s')
% ylabel('rad/s^2')
% axis([0 time -25 25]);
% title('angular acceleration of fishtail relative to fish head');
% 
% figure (14)
% plot(x_axis,X_0_2o_t); %�����н����ٶ�
% xlabel('time/s')
% ylabel('m/s^2')
% axis([0 time -0.2 0.2]);
% title('traveling acceleration of fish');
% 
% figure (15)
% plot(x_axis,theta_10_g_t);%�����������������������ͷ���Ƕ�
% xlabel('time/s')
% ylabel('rad')
% axis([0 time -0.2 0.2]);
% title('given angle of fish body relative to fish head');
% 
% figure (16)
% plot(x_axis,theta_21_g_t);%����������β���������������Ƕ�
% xlabel('time/s')
% ylabel('rad')
% axis([0 time -1 1]);
% title('given angle of fishtail relative to fish body');
% 
% figure (17)
% plot(x_axis,theta_20_g_t);%����������β�����������ͷ���Ƕ�
% xlabel('time/s')
% ylabel('rad')
% axis([0 time -1 1]);
% title('given angle of fishtail relative to fish body');
% 
% figure (18)
% plot(x_axis,M_1_t);%��һ���ؽڵ�ת��
% xlabel('time/s')
% ylabel('N*m')
% axis([0 time -0.4 0.4]);
% title('torque of the first joint');
% 
% figure (19)
% plot(x_axis,M_2_t);%�ڶ����ؽڵ�ת��
% xlabel('time/s')
% ylabel('N*m')
% axis([0 time -0.2 0.2]);
% title('torque of the second join');
% 
% figure (20)
% plot(x_axis,F_1_s_t);%�������岿��������������Ϊ��
% xlabel('time/s')
% ylabel('N')
% axis([0 time -10 25]);
% title('fish body force, scalar, positive');
% 
% figure (21)
% plot(x_axis,F_2_s_t);%����β��������������Ϊ��
% xlabel('time/s')
% ylabel('N')
% axis([0 time -10 25]);
% title('fishtail force, scalar, positive');
% 
% figure (22)
% plot(x_axis,M_1_t,x_axis,theta_10_1o_t);%%,'r',x_axis,x_zeros); %��һ�ؽ�ת�ؼ������岿���������ͷ���Ľ��ٶ�
% legend('M_1','\theta_{10}dot')
% xlabel('time/s')
% ylabel('M/N*m w/rad^-1')
% axis([0 time -1.2 1.2]);
% title('torque of the first joint and angular velocity of fish body relative to fish head');
% 
% 
% figure (23)
% plot(x_axis,M_2_t,x_axis,theta_20_1o_t);%%,'r',x_axis,x_zeros); %�ڶ��ؽ�ת�ؼ���β�������������Ľ��ٶ�
% legend('M_2','\theta_20dot')
% xlabel('time/s')
% ylabel('M/N*m w/rad*s^-1')
% axis([0 time -4 4]);
% title('torque of the second joint and angular velocity of fishtail relative to fish boay');
% 
% 
% 
% figure (24)
% plot(x_axis,M_1_t,x_axis,theta_10_t);%%,'r',x_axis,x_zeros); %��һ�ؽ�ת�ؼ������岿���������ͷ���ĽǶ�
% legend('M_1','\theta_{10}')
% xlabel('time/s')
% ylabel('M/N*m ��/rad*s^-1')
% axis([0 time -0.3 0.3]);

% 
% 
% 
% figure (25)
% plot(x_axis,M_2_t,x_axis,theta_21_t);%%,'r',x_axis,x_zeros);  %�ڶ��ؽ�ת�ؼ���β������������岿�ֵĽǶ�
% legend('M_2','\theta_{21}')
% xlabel('time/s')
% ylabel('M/N*m ��/rad*s^-1')
% axis([0 time -0.7 0.7]);
% % title('torque of the first joint and angule of fishtail relative to fish body');



