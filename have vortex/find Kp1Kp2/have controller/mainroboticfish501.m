%%%主程序，调用机器鱼的子程序来寻找最合适的参数f,kp1,kp2
%%%寻找的标准为使效率最高的那一组参数，从而可以用来画出拟合的曲线。
i_x=0;
j_x=0;
k_x=0;
f_x=0;
H_1_x=0;
H_2_x=0;
H_1_x_t=zeros(1,600);
H_2_x_t=zeros(1,600);
f_x_t=zeros(1,4);
f_x_t(1)=1;
f_x_t(2)=1.5;
f_x_t(3)=1.8;
f_x_t(4)=2.0;
outyita_x_t=zeros(60,600);

    for H_1_x=10:1:15
        j_x=j_x+1;
        i_x=0;
        for H_2_x=10:1:15
            i_x=i_x+1;
           outyita_x_t(j_x,i_x)=roboticfish501(2,7,23,H_1_x,H_2_x);
         
           H_1_x_t(i_x)=H_1_x;
           H_2_x_t(i_x)=H_2_x;
             disp(['f_x = ',num2str(f_x)]);
            disp(['Kp_1_x = ',num2str(H_1_x)]);
            disp(['Kp_2_x = ',num2str(H_2_x)]);
        end;
    end;

