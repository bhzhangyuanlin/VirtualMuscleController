%%%�����򣬵��û�������ӳ�����Ѱ������ʵĲ���f,kp1,kp2
%%%Ѱ�ҵı�׼ΪʹЧ����ߵ���һ��������Ӷ���������������ϵ����ߡ�
i_x=0;
j_x=1;
k_x=0;
f_x=0;
Kp_1_x=0;
Kp_2_x=0;
Kp_1_x_t=zeros(1,40);
Kp_2_x_t=zeros(1,40);
f_x_t=zeros(1,1000);
outyita_x_t=zeros(600,40);


        for Kp_1_x=1:1:30
            i_x=i_x+1;
           outyita_x_t(j_x,i_x)=roboticfish209(Kp_1_x,16);
            Kp_1_x_t(i_x)=Kp_1_x;
            Kp_2_x_t(i_x)=Kp_2_x;
            disp(['Kp_1_x = ',num2str(Kp_1_x)]);
            disp(['Kp_2_x = ',num2str(Kp_2_x)]);
        end;
