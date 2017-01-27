[p_old,theta,differ,line_out,name] = plot_phase_old('\1125\hato_static',200);
[image_height,image_width] = size(theta);
unwrapped_phase_temp = zeros(image_height,image_width);
medfilt = [5 5];


folder_name = 'C:\Users\t2ladmin\Documents\MATLAB\phase shift\phase_image';
file_name = [folder_name, '\', name];
cd(file_name);

%%%%%%%%%%%%%%%%%%%%%%微分情報を格納,外れ値を除去
temp_p1= p_old(:,1);
A = find(temp_p1<=0);
temp_p2= p_old(:,2);
B = find(temp_p2>=0);
temp_p = p_old;
for i = 1:length(A)
    temp_p(A(i),:)=[0,0];
end
for i = 1:length(B)
    temp_p(B(i),:)=[0,0];
end
temp_p(temp_p==0)=[];
temp_p = reshape(temp_p,[length(temp_p)/2 2]);
%%%%%%%%%%%%%%%%%%%%%%%%絶対位相が0となるピクセルを推定
phase_0 = round(-mean(temp_p(:,2))/mean(temp_p(:,1)))*ones(image_height,1);
p = horzcat(mean(temp_p(:,1))*ones(image_height,1),mean(temp_p(:,2))*ones(image_height,1));
%phase_0 = 132*ones(image_height,1);
%p = horzcat(0.0007*ones(image_height,1),-0.0926*ones(image_height,1));

%%%%%%%%%%%%%%%%%位相の微分情報をプロット
%グラフ4
figure;
plot(1:image_width-1,differ(line_out,:));
hold on
plot(1:image_width-1,polyval(p(line_out,:),1:image_width-1),'r'); %位相の微分の線形近似をプロット
title(['位相の微分',num2str(line_out)],'FontSize',16);
xlim([0,image_width]);
ylim([-0.5,0.5]);
set(gca,'FontSize',16);
saveas(gcf,['derivative_',num2str(line_out),'.fig'])


%%%%%%%%%%%%%%%%%微分情報から位相を表す2次関数を導出       
for i = 1:image_height
    if p(i,1) > 0 
        %%位相接続
        for j = 1:image_width
            unwrapped_phase_temp(i,j) = 0.5*p(i,1)*(phase_0(i)-j)^2;
            %{
            if 0<phase_0(i) && phase_0(i)<image_width+1
                unwrapped_phase_temp(i,j) = unwrapped_phase_temp(i,j);
            end
            %}
        end
    else
        unwrapped_phase_temp(i,:) = NaN(1,image_width);
    end
end
          
%%%%%%%%%%%%%%%%%%%%%%%%%位相接続%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_map = unwrapped_phase_temp-theta;
temp_out = temp_map/(2*pi);
index = round(temp_out); %位相接続の際のインデックスを格納
unwrapped_phase = index*2*pi + theta; %絶対位相を格納
%プロジェクタとの対応をとる
pixel = 200*sqrt(abs(unwrapped_phase));

%%%%%%%%%%%%%%%%%%%%%二次関数をプロット
%%グラフ5
figure;
plot(1:image_width,theta(line_out,:));
hold on;
plot(1:image_width,unwrapped_phase_temp(line_out,:),'r');
title(['位相',num2str(line_out)],'FontSize',16);
xlim([0,image_width]);
set(gca,'FontSize',16);
saveas(gcf,['wrap_and_unwrap_',num2str(line_out),'.fig'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%絶対位相をプロット
%%グラフ6
figure;
pixel_lim = [0 1024];
%pixel = medfilt2(pixel, medfilt);
imagesc(pixel,pixel_lim);
colorbar;
hold on;
plot(1:image_width,line_out,'r');
title('絶対位相','Fontsize',16);
set(gca,'FontSize',16);
saveas(gcf,'unwrappedphase_1.fig')

%%%%%%%%%%%%%%%%%%%%%位相接続時の基準(二次関数)とラップされた位相の差をプロット
%%グラフ7
figure;
plot(1:image_width,temp_out(line_out,:),'r');
hold on;
plot(1:image_width,index(line_out,:),'g');
hold on;
plot(1:image_width,theta(line_out,:));
xlim([0,image_width]);
title(['index',num2str(line_out)],'Fontsize',16);
saveas(gcf,['index_',num2str(line_out),'.fig'])

%%%%%%従来の位相接続
index2 = zeros(image_height,image_width);
for i = 1:image_height
    for j = 1:image_width-1
        if differ(i,j)< -pi
            index2(i,j+1:image_width) = index2(i,j+1:image_width) + ones(1,image_width-j);
        elseif isnan(theta(i,j))==1
            index2(i,j)=NaN;
        end
    end
end
unwrapped_phase2 = theta + index2*2*pi;
pixel2 = 200*sqrt(abs(unwrapped_phase2));
%%%%%%従来の位相接続手法で求めた絶対位相をプロット
%%%グラフ8
figure;
imagesc(pixel2,pixel_lim);
colorbar;
hold on;
plot(1:image_width,line_out,'r');
title('絶対位相（従来法）','Fontsize',16);
set(gca,'FontSize',16);
saveas(gcf,'unwrappedphase_2.fig')

figure
plot(1:image_width,pixel(line_out,:));
title('プロジェクタピクセル','Fontsize',16);
xlabel('カメラ列','Fontsize',16)
ylabel('ピクセル','Fontsize',16)
ylim([0 1024]);
set(gca,'FontSize',16);
saveas(gcf,['pixel_info',num2str(line_out),'.fig'])



cd('C:\Users\t2ladmin\Documents\MATLAB\phase shift');