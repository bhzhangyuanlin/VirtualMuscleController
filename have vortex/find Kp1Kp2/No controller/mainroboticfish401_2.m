%%%主程序，调用机器鱼的子程序来寻找最合适的参数f,kp1,kp2
%%%寻找的标准为使效率最高的那一组参数，从而可以用来画出拟合的曲线。
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