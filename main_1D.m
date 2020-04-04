    clc;
    clear all;
    clear global;
    close all;
    
    oned_testmode = 'compara_esprit';
    
    switch oned_testmode
        case 'compara_esprit'
            theta0 = [-60 -45 -30 -15 0 20 40 60];
            [rmse1]=conjugate_esprit(theta0,8);
            [rmse2]=espirit(theta0,8);
            figure('Color','white');
            plot(-10:1:9,rmse1,'g*-',-10:1:9,rmse2,'ro-');
            grid on;
            xlabel('SNR/dB');
            ylabel('RMSE/degree');
            title('共轭TLS-ESPRIT与TLS-ESPRIT算法在不同信噪比下的均方根误差');
            legend('conjugate esprit','esprit');
        case 'conjugate_esprit'
            theta0 = [-60 -45 -30 -15 0 20 40  60];
            [rmse]=conjugate_esprit(theta0,8);
            figure('Color','white');
            plot(rmse,'*-');
            grid on;
            xlabel('SNR/dB');
            ylabel('RMSE/degree');
            title('共轭ESPRIT算法在不同信噪比下的均方根误差');
        case 'espirit'
            theta0 = [-60 -45 -30  0 20 40  60];
            [rmse]=espirit(theta0,8);
            figure('Color','white');
            plot(rmse,'*-');
            grid on;
            xlabel('SNR/dB');
            ylabel('RMSE/degree');
            title('ESPRIT算法在不同信噪比下的均方根误差');
        case 'music_and_espirit'
            theta0 = [10 30 60];
            theta=linspace(-90,90,361);
            [abs_p,peak_ang,snr1,rmse1] = music(theta0,8,0.5,'not_multi');
            [rmse2] = espirit(theta0,8);
            figure('Color','white');
            plot(snr1,rmse1,'r*-',snr1,rmse2,'b*-');
            title('MUSIC算法与ESPRIT算法不同信噪比下测角均方误差');
            xlabel('SNR/dB');ylabel('RMSE/degree');
            legend('MUSIC','ESPIRIT');
        case 'rootmusic_and_espirit'
            theta0 = [0 30 60];
            [snr1,rmse1] = root_music(theta0,8,0.5);  
            [rmse2] = espirit(theta0,8);
            figure('Color','white');
            plot(snr1,rmse1,'ro-',snr1,rmse2,'g*-');
            grid on;
            title('ROOT-MUSIC算法与ESPRIT算法不同信噪比下测角均方误差');
            xlabel('SNR/dB');ylabel('RMSE/degree');
            legend('ROOT-MUSIC','ESPIRIT');
        case 'nor_and_ss_music'
            theta0 = [10 30 60];
            theta=linspace(-90,90,361);
            [abs_p1,peak_ang1] = music(theta0,8,0.5,'multi_path');
            abs_p2 = ss_music(theta0,8,0.5,'multi_path');
            figure('Color','white');
            h=plot(theta,abs_p1,'g',theta,abs_p2,'r');
            set(h,'Linewidth',2)
            xlabel('angle (degree)')
            ylabel('abs of P (dB)')
            axis([-90 90 -60 0]);
            set(gca, 'XTick',[-90:30:90]);
            grid on;
            legend('MUSIC','空间平滑MUSIC');
            title({['MUSIC算法和空间平滑MUSIC算法在多径传输下的性能比较']...
                    ['来波方向为:' num2str(theta0) '信源相同的方向为:' num2str(theta0(1:2))]});
        case 'root_music'
            theta0 = [45 30 60];
            [snr,rmse] = root_music(theta0,8,0.5);  
            figure('Color','white');
            plot(snr,rmse,'ro-');
            title('ROOT-MUSIC算法不同信噪比下测角均方误差');
            xlabel('SNR/dB');ylabel('RMSE/degree');
        case 'ss_music'
            theta=linspace(-90,90,361);
            theta0 = [0 30 60];
            abs_p = ss_music(theta0,8,0.5,'multi_path');
            figure('Color','white');
            h=plot(theta,abs_p);
            set(h,'Linewidth',2);
            xlabel('angle (degree)');
            ylabel('abs of P (dB)');
            axis([-90 90 -60 0]);
            set(gca, 'XTick',[-90:30:90], 'YTick',[-60:10:0]);
            grid on;
            title('平滑MUSIC功率谱');
        case 'music'
            theta0 = [10 20 30 40];
            [abs_p,peak_ang,snr,rmse] = music(theta0,8,0.5,'not_multi');
            theta=linspace(-90,90,361);
            figure('Color','white');
            h=plot(theta,abs_p);
            hold on;
            plot(peak_ang,abs_p(peak_ang*2+181),'ro');
            set(h,'Linewidth',2);
            xlabel('angle (degree)');
            ylabel('abs of P (dB)');
            axis([-90 90 -60 0]);
            set(gca, 'XTick',[-90:30:90]);
            grid on;
            title( {['MUSIC功率谱,来波方向为:' num2str(theta0)]...
                ['谱峰搜索结果为:' num2str(peak_ang)]});
            figure('Color','white');
            plot(snr,rmse,'ro-');
            title('MUSIC算法不同信噪比下测角均方误差');
            xlabel('SNR');ylabel('RMSE');
        case 'normal_theta_and_direction'
            %% 普通波束形成的波束宽度随方向的变化
            theta=linspace(-90,90,200);
            d_lamda = 1/2;
            theta0 = [0 30 60 80];
            [abs_f1,abs_p1] = normal(theta0(1),16,d_lamda);
            [abs_f2,abs_p2] = normal(theta0(2),16,d_lamda);
            [abs_f3,abs_p3] = normal(theta0(3),16,d_lamda);
            [abs_f4,abs_p4] = normal(theta0(4),16,d_lamda);

            figure('Color','white');
            subplot(411);
            plot(theta,abs_f1,'r');grid on;
            title(['16阵元线阵方向图，指向方向为',num2str(theta0(1)) '度']);
            xlabel('theta/degree');ylabel('abs of F');
            subplot(412);
            plot(theta,abs_f2,'r');grid on;
            title(['16阵元线阵方向图，指向方向为',num2str(theta0(2)) '度']);
            xlabel('theta/degree');ylabel('abs of F');
            subplot(413);
            plot(theta,abs_f3,'r');grid on;
            title(['16阵元线阵方向图，指向方向为',num2str(theta0(3)) '度']);
            xlabel('theta/degree');ylabel('abs of F');
            subplot(414);
            plot(theta,abs_f4,'r');grid on;
            title(['16阵元线阵方向图，指向方向为',num2str(theta0(4)) '度']);
            xlabel('theta/degree');ylabel('abs of ');
        case 'normal_resolution'
            %% 普通波束形成的分辨力问题 
            theta=linspace(-90,90,200);
            d_lamda = 1/2;
            theta1 = [0 30 60];
            [abs_f1,abs_p1] = normal(theta1,8,d_lamda);
            theta2 = [0 10 60];
            [abs_f2,abs_p2] = normal(theta2,8,d_lamda);
            theta3 = [0 50 60];
            [abs_f3,abs_p3] = normal(theta3,8,d_lamda);
            figure('Color','white');
            subplot(311);
            plot(theta,abs_p1,'r');grid on;
            title(['8阵元线阵功率谱，来波方向为',num2str(theta1) '度']);
            xlabel('theta/degree');ylabel('abs of P');
            subplot(312);
            plot(theta,abs_p2,'r');grid on;
            title(['8阵元线阵功率谱，来波方向为',num2str(theta2) '度']);
            xlabel('theta/degree');ylabel('abs of P');
            subplot(313);
            plot(theta,abs_p3,'r');grid on;
            title(['8阵元线阵功率谱，来波方向为',num2str(theta3) '度']);
            xlabel('theta/degree');ylabel('abs of P');

        case 'theta_and_array_numb'
            %% 测试阵元数目对波束宽度的影响
            theta0 = 0;
            d_lamda = 1/2;
            [abs_f1,abs_p1] = normal(theta0,8,d_lamda);
            [abs_f2,abs_p2] = normal(theta0,16,d_lamda);
            [abs_f3,abs_p3] = normal(theta0,32,d_lamda);
            theta=linspace(-90,90,200);
            figure('Color','white');
            subplot(211);
            plot(theta,abs_f1,'g');hold on;
            plot(theta,abs_f2,'r');hold on;
            plot(theta,abs_f3,'b');hold on;
            grid on;
            xlabel('theta/degree')
            ylabel('abs of F')
            title([ '均匀线阵方向图','指向' num2str(theta0) '度' ]);
            legend('8阵元','16阵元','32阵元');
            subplot(212);
            plot(theta,abs_p1,'g');hold on;
            plot(theta,abs_p2,'r');hold on
            plot(theta,abs_p3,'b');hold on;
            grid on;
            xlabel('theta/degree')
            ylabel('abs of P')
            title([ '均匀线阵功率谱','来波方向为' num2str(theta0) '度']);
            legend('8阵元','16阵元','32阵元');

        case 'theta_and_array_span'
            %% 测试阵元间距对波束宽度的影响
            theta0 = 0;
            d_lamda = 1/2;
            [abs_f1,abs_p1] = normal(theta0,8,d_lamda/2);
            [abs_f2,abs_p2] = normal(theta0,8,d_lamda);
            [abs_f3,abs_p3] = normal(theta0,8,d_lamda/2*3);
            [abs_f4,abs_p4] = normal(theta0,8,d_lamda*2*1.5);
            theta=linspace(-90,90,200);
            figure('Color','white');
            subplot(411);
            plot(theta,abs_f1,'g');
            grid on;
            xlabel('theta/degree')
            ylabel('abs of F')
            title([ '阵元间距' num2str(d_lamda/2) 'lamda均匀线阵方向图','指向' num2str(theta0) '度' ]);
            subplot(412);
            plot(theta,abs_f2,'r');
            grid on;
            xlabel('theta/degree')
            ylabel('abs of F')
            title([ '阵元间距' num2str(d_lamda) 'lamda均匀线阵方向图','指向' num2str(theta0) '度' ]);
            subplot(413);
            plot(theta,abs_f3,'b');
            grid on;
            xlabel('theta/degree')
            ylabel('abs of F')
            title([ '阵元间距' num2str(d_lamda/2*3) 'lamda均匀线阵方向图','指向' num2str(theta0) '度' ]);
            subplot(414);
            plot(theta,abs_f4,'k');
            grid on;
            xlabel('theta/degree')
            ylabel('abs of F')
            title([ '阵元间距' num2str(d_lamda*2*1.5) 'lamda均匀线阵方向图','指向' num2str(theta0) '度' ]);

        case 'normal_and_capon'
            %% 对比一维普通波束形成与capon
            theta0 = [30,60];
            d_lamda = 1/2;
            [abs_f1,abs_p1] = normal(theta0,8,d_lamda);
            [abs_f2,abs_p2] = capon(theta0,16,d_lamda);
            theta=linspace(-90,90,200);
            figure('Color','white');
            subplot(211);
            plot(theta,abs_p1,'r');hold on;
            grid on;
            xlabel('theta/degree');
            ylabel('abs of P');
            title([ '均匀线阵功率谱(普通波束形成)','来波方向为' num2str(theta0) '度']);
            subplot(212);
            plot(theta,abs_p2,'r');hold on;
            grid on;
            xlabel('theta/degree');
            ylabel('abs of P');
            title([ '均匀线阵功率谱(Capon)','来波方向为' num2str(theta0) '度']);
            theta0 = [30,40];
            d_lamda = 1/2;
            [abs_f1,abs_p1] = normal(theta0,8,d_lamda);
            [abs_f2,abs_p2] = capon(theta0,16,d_lamda);
            figure('Color','white');
            subplot(211);
            plot(theta,abs_p1,'r');hold on;
            grid on;
            xlabel('theta/degree')
            ylabel('abs of P')
            title([ '均匀线阵功率谱(普通波束形成)','来波方向为' num2str(theta0) '度']);
            subplot(212);
            plot(theta,abs_p2,'r');hold on;
            grid on;
            xlabel('theta/degree')
            ylabel('abs of P')
            title([ '均匀线阵功率谱(Capon)','来波方向为' num2str(theta0) '度']);
        otherwise
            disp('Please input mode!');
    end