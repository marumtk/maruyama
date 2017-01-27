%取得したい画像があるフォルダ名を入力する
%位相の回帰係数，各ピクセルの位相，x方向の位相の微分を出力
function [p,theta,differ_raw,line_out,name] = plot_phase_old(name,line_out)
    close all;
    %%%%%%%%%%
    %%%%%%%%%%
    
    folder_name = 'C:\Users\t2ladmin\Documents\MATLAB\phase shift\phase_image';
    file_name = [folder_name, '\', name];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%フォルダ内に移動し，画像のみを読み込み
    cd(file_name);
    imgFileList = dir('*.png');  %%読み込む画像はtifかpngなのでそのつど書き換えが必要
    imgFileList.name
    
    num = length(imgFileList); %画像の数
    image_width = 512; %カメラのピクセル数(x方向)
    image_height = 512; %カメラのピクセル数(y方向)
    pixel = image_width * image_height; %カメラの全ピクセル数
    Luminance = zeros(pixel,num); %全画像，全ピクセルの輝度を格納する行列

    %実際に輝度を記録する
    for i = 1:num
        img_name = char(imgFileList(i).name);
        temp = imread(img_name);  %カメラの撮影画像の読み込み
        %カラー画像の場合2値化する
        if(numel(temp) == pixel*3)
            temp = rgb2gray(temp);
        end
        temp = temp(:);
        Luminance(:,i) = temp;
        clear temp
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%輝度の記録終了
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%位相の計算
    Luminance = Luminance';
    M = [1,-1/2,-sqrt(3)/2;1,1,0;1,-1/2,+sqrt(3)/2];
    temp = M\Luminance;
    clear M;
    I_mod = sqrt(temp(2,:).^2+temp(3,:).^2); %振幅の導出
   
    theta = zeros(1,pixel); %位相を格納
    %位相を計算していく
    for i = 1:pixel
        if I_mod(i)<10
            theta(i)=-10; %振幅が小さいとき，正確な位相の計算は不可能
        else
            sin = temp(2,i)./I_mod(i);
            cos = temp(3,i)./I_mod(i);
            if abs(cos)<=1 && abs(sin)<=1
                theta(i) = atan2(sin,cos);
            else
                theta(i)=-10;
            end
        end
    end
    
    theta = reshape(theta,[image_height image_width]); %-piからpiの範囲で位相を求める
    medfilt = [5 5];
    theta = medfilt2(theta, medfilt);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%位相の計算終了
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%位相の微分に着目
    differ = zeros(image_height,image_width-1); %thetaの差（x方向の微分）を格納
    differ_raw = zeros(image_height,image_width-1); %thetaの差（x方向の微分）の生データを格納
    p = zeros(image_height,2); %thetaの微分の回帰を行ごとに格納
    cd('C:\Users\t2ladmin\Documents\MATLAB\phase shift');
    for line = 1:image_height
        %位相の空間方向の微分を計算        
        differ_raw(line,:) = diff(theta(line,:));
        differ(line,:) = sgolayfilt(differ_raw(line,:),1,17);
    
        %外れ値を除外し，位相の微分を線形近似
        differ_temp = zeros(2,image_width-1);
        num = 1;
        for i = 1:image_width-1
            if differ(line,i)>0 && differ(line,i)<0.3 
                differ_temp(1,num) = i;
                differ_temp(2,num) = differ(line,i);
                num = num+1;
            end
        end
        p(line,:) = polyfit(differ_temp(1,1:num-1),differ_temp(2,1:num-1),1);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%位相の微分の計算と回帰終了
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ある行についての位相や輝度をプロットする
    %%%%%%位相のプロット
    %グラフ1
    cd(file_name);
    figure;
    plot(1:image_width,theta(line_out,:));
    title(['位相',num2str(line_out)],'FontSize',16);
    xlim([0,image_width]);
    ylim([-pi,pi]);
    set(gca,'FontSize',16);
    saveas(gcf,['phase_',num2str(line_out),'.fig'])
    
    %%%%%%位相の空間方向の微分をプロット
    %{
    figure;
    plot(1:image_width-1,differ(line_out,:));
    hold on
    plot(1:image_width-1,polyval(p(line_out,:),1:image_width-1),'r'); %位相の微分の線形近似をプロット
    title(['位相の微分',num2str(line_out)],'FontSize',16);
    xlim([0,image_width]);
    ylim([-0.5,0.5]);
    set(gca,'FontSize',16);
    saveas(gcf,['derivative_',num2str(line_out),'.fig'])
    %}
    
    %%%%%%輝度をプロット，sin波になってればOK
    %グラフ2
    figure;
    sin_wave = reshape(Luminance(2,:), [image_height image_width]);
    plot(1:image_width,sin_wave(line_out,:));
    title(['輝度',num2str(line_out)],'FontSize',16);
    xlim([0,image_width]);
    ylim([0,256]);
    set(gca,'FontSize',16);
    saveas(gcf,['I_',num2str(line_out),'.fig'])
  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%実際の位相復元画像を表示．不要な場合は以下をコメントアウト
    %グラフ3
    figure;
    theta_lim = [-pi pi];
    imagesc(theta,theta_lim);
    colorbar;
    hold on;
    %plot(1:image_width,line_out,'r');
    title('位相マップ','FontSize',16);
    set(gca,'FontSize',16);
    saveas(gcf,'phasemap.fig')
   
    figure;
    plot(1:image_width-1,differ(line_out,:));
    hold on
    plot(1:image_width-1,differ_raw(line_out,:),'r'); %位相の微分の線形近似をプロット
    title(['位相の微分',num2str(line_out)],'FontSize',16);
    xlim([0,image_width]);
    ylim([-0.5,0.5]);
    set(gca,'FontSize',16);

    %%元のディレクトリに戻る
    cd('C:\Users\t2ladmin\Documents\MATLAB\phase shift');
end