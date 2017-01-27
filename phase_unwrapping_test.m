[p_old,theta,differ,line_out,name] = plot_phase_old('\1220\hand_proposed',350);
[image_height,image_width] = size(theta);
unwrapped_phase_temp = zeros(image_height,image_width);
phase_0 = zeros(image_height,1); %位相の微分が0となる点を推定
p = zeros(image_height,2); %平均微分情報を格納
medfilt = [10 10];


folder_name = 'C:\Users\t2ladmin\Documents\MATLAB\phase shift\phase_image';
file_name = [folder_name, '\', name];
cd(file_name);

%%%%%%%%%%%%%%%%%%%%%%微分情報を格納,外れ値を除去
temp_p1= p_old(:,1);
A = find(temp_p1<=0);
temp_p = p_old;
Ind = 1:768; %%生成した1024*768の画像を用いる時はIndも中身も変わる

%%生の画像データを用いる時は以下をコメントアウト
%外れ値除去
for i = 1:length(A)
    temp_p(A(i),:)=[NaN,NaN];
end

%%ピクセル数が少ない行を除去
for i = 1:512
    Temp = find(theta(i,:)~=-10);
    if length(Temp) < 50
        temp_p(i,:) = [NaN,NaN];
    end
end
Ind = find(isnan(temp_p(:,1))==0); %値が入っている行を格納
temp_p(isnan(temp_p)==1)=[];
temp_p = reshape(temp_p,[length(temp_p)/2 2]);

%端を除去
exception = ceil(length(temp_p)*0.1);
out = length(temp_p);
for i = 1:exception
    temp_p(i,:)=[0,0];
    temp_p(out-i+1,:)=[0,0];
    Ind(i,:)=0;
    Ind(out-i+1,:)=0;
end
temp_p(temp_p==0)=[];
temp_p = reshape(temp_p,[length(temp_p)/2 2]);
Ind(Ind==0)=[];

result_len = length(temp_p(:,1));
temp_p(:,1) = medfilt1(temp_p(:,1),30);
temp_p(:,2) = medfilt1(temp_p(:,2),30);
temp_p(:,1) = mean(temp_p(:,1))*ones(result_len,1);
temp_p(:,2) = mean(temp_p(:,2))*ones(result_len,1);
phase_0_temp = -temp_p(:,2)./temp_p(:,1);
%%生の画像データを用いる時はここまでをコメントアウト

%%%%%%%%%%%%%%%%%%%%%%%%絶対位相が0となるピクセルを推定
Ind = vertcat(1,Ind);
Ind = vertcat(Ind,image_width);
for i = 1:length(Ind)-1
    if i~=length(Ind)-1
        for j = Ind(i):Ind(i+1)
            phase_0(j) = phase_0_temp(i);
            p(j,1) = temp_p(i,1);
            p(j,2) =- p(j,1)*phase_0(j);
        end
    else
        for j = Ind(i):Ind(i+1)
            phase_0(j) = phase_0_temp(length(temp_p));
            p(j,1) = temp_p(length(temp_p));
            p(j,2) =- p(length(temp_p),1)*phase_0(length(temp_p));
        end
    end    
end    
%phase_0 = round(-mean(temp_p(:,2))/mean(temp_p(:,1)))*ones(image_height,1);
%p = horzcat(mean(temp_p(:,1))*ones(image_height,1),mean(temp_p(:,2))*ones(image_height,1));
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
            if theta(i,j) ~= -10
                unwrapped_phase_temp(i,j) = 0.5*p(i,1)*(phase_0(i)-j)^2;
            else
                unwrapped_phase_temp(i,j) = -10;
            end
        end
    else
        unwrapped_phase_temp(i,:) = -10*ones(1,image_width);
    end
end
          
%%%%%%%%%%%%%%%%%%%%%%%%%位相接続%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_map = unwrapped_phase_temp-theta;
temp_out = temp_map/(2*pi);
index = round(temp_out); %位相接続の際のインデックスを格納
unwrapped_phase = index*2*pi + theta; %絶対位相を格納
unwrapped_phase = medfilt2(unwrapped_phase,medfilt);
unwrapped_phase(unwrapped_phase==-10) = NaN;
%プロジェクタとの対応をとる
pixel = 200*sqrt(abs(unwrapped_phase));
%誤差が大きいときは二次関数を利用
%pixel_temp = 200*sqrt(abs(unwrapped_phase_temp));
%pixel(abs(pixel_temp-pixel)>1) = pixel_temp(abs(pixel_temp-pixel)>1);

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
colormap(jet(256));
axis image;
hold on;
%plot(1:image_width,line_out,'r');
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

%%%%%%%%%%%%%%%%%%%%%プロジェクタのピクセルをプロット
%%グラフ8
figure
plot(1:image_width,pixel(line_out,:));
title('プロジェクタピクセル','Fontsize',16);
xlabel('カメラ列','Fontsize',16)
ylabel('ピクセル','Fontsize',16)
xlim([0 image_width]);
set(gca,'FontSize',16);
saveas(gcf,['pixel_info',num2str(line_out),'.fig'])


cd('C:\Users\t2ladmin\Documents\MATLAB\phase shift');