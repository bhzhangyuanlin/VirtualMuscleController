%%%�����򣬵��û�������ӳ�����Ѱ������ʵĲ���f,kp1,kp2
%%%Ѱ�ҵı�׼ΪʹЧ����ߵ���һ��������Ӷ���������������ϵ����ߡ�
i_x=0;
j_x=1;
k_x=0;
f_x=0;
H_1_x=0;
H_2_x=9;
H_1_x_t=zeros(1,40);
H_2_x_t=zeros(1,40);
f_x_t=zeros(1,1000);
outyita_x_t=zeros(600,40);


        for H_1_x=-10:1:0
            i_x=i_x+1;
           outyita_x_t(j_x,i_x)=roboticfish212(H_1_x,9);
            H_1_x_t(i_x)=H_1_x;
            H_2_x_t(i_x)=H_2_x;
            disp(['H_1_x = ',num2str(H_1_x)]);
            disp(['H_2_x = ',num2str(H_2_x)]);
        end;
