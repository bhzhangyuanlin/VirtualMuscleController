%%%�����򣬵��û�������ӳ�����Ѱ������ʵĲ���f,kp1,kp2
%%%Ѱ�ҵı�׼ΪʹЧ����ߵ���һ��������Ӷ���������������ϵ����ߡ�
i_x=0;
j_x=0;
k_x=0;
f_x=0;
Kp_1_x=0;
Kp_2_x=0;
Kp_1_x_t=zeros(1,40);
Kp_2_x_t=zeros(1,40);
f_x_t=zeros(1,1000);
outyita_x_t=zeros(600,40);
for f_x=1.5:0.3:1.8
    for Kp_1_x=1:1:10
        j_x=j_x+1;
        i_x=0;
        for Kp_2_x=1:1:17
            i_x=i_x+1;
           outyita_x_t(j_x,i_x)=roboticfish401(f_x,Kp_1_x,Kp_2_x);
         
            Kp_1_x_t(i_x)=Kp_1_x;
            Kp_2_x_t(i_x)=Kp_2_x;
             disp(['f_x = ',num2str(f_x)]);
            disp(['Kp_1_x = ',num2str(Kp_1_x)]);
            disp(['Kp_2_x = ',num2str(Kp_2_x)]);
        end;
    end;
end;