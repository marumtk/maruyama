[p,theta,differ,line_out,name,gap] = plot_phase('\1125\hand_static_noref',50);
[image_height,image_width] = size(theta); %カメラピクセルのサイズ

unwrapped_phase_temp = zeros(image_height,image_width); %仮の絶対位相を格納(二次関数近似の結果)
partition = 1; %微分情報を求める際の分割数(ブロック数)
block_height = image_height/partition; %微分情報を求める際の一ブロックのサイズ
partition_noise = zeros(partition,2); %一ブロックあたりのノイズ数
medfilt = [5 5]; %メディアンフィルタ

folder_name = 'C:\Users\t2ladmin\Documents\MATLAB\phase\phase_image';
file_name = [folder_name, '\', name];
cd(file_name);


%%%%%%%%%%%%%%%%%%%%%%微分情報を格納,外れ値を除去
temp_p1= p(:,1);
A = find(temp_p1<=0);
A1 = find(temp_p1>1);
temp_p2= p(:,2);
B = find(temp_p2>=0);
B1 = find(temp_p2<-1);
temp_p = p;
for i = 1:length(A)
    temp_p(A(i),:)=[0,0];
end
for i = 1:length(A1)
    temp_p(A1(i),:)=[0,0];
end
for i = 1:length(B)
    temp_p(B(i),:)=[0,0];
end
for i = 1:length(B1)
    temp_p(B1(i),:)=[0,0];
end


%%%%%%%%%%%%%%%%%%%%%二次関数の頂点となるピクセルを列ごとに推定
phase_0 = zeros(image_height,1); %theta^2を表す二次関数の頂点となるピクセルを格納

%%%%%%%%%%%%基準平面を用意して平行ステレオとなるように撮像するならこのアルゴリズムが使える
%{
for i = 1:image_height
    temporal = find(theta(i,:)~=0);
    if isempty(temporal)==1
        phase_0(i,1)=0;
    else
        phase_0(i,1)=temporal(1);
    end
end
%}

%%%%%%%%%%%%%%%%%ノイズの数をブロックごとに算出
for i = 1:partition
    num1 = find(temp_p(1+(i-1)*block_height:i*block_height,1)==0);
    num2 = find(gap(1+(i-1)*block_height:i*block_height,1)==0);
    partition_noise(i,1) = length(num1); %微分情報のノイズの数
    partition_noise(i,2) = length(num2); %二次関数の頂点についてのノイズの数
end

%%%%%%%%%%%%%%%%%微分情報の平均値をブロックごとに導出        
for i = 1:partition
    p(1+(i-1)*block_height:i*block_height,1) = mean(temp_p(1+(i-1)*block_height:i*block_height,1))*(block_height/(block_height-partition_noise(i,1)))*ones(block_height,1);
    p(1+(i-1)*block_height:i*block_height,2) = mean(temp_p(1+(i-1)*block_height:i*block_height,2))*(block_height/(block_height-partition_noise(i,1)))*ones(block_height,1);
    phase_0(1+(i-1)*block_height:i*block_height,1) = round(-p(i*block_height,2)/p(i*block_height,1))*ones(block_height,1);
    if partition_noise(i,2) ~= block_height
        gap(1+(i-1)*block_height:i*block_height,1) = mean(gap(1+(i-1)*block_height:i*block_height,1))*(block_height/(block_height-partition_noise(i,2)))*ones(block_height,1);
    else
        gap(1+(i-1)*block_height:i*block_height,1) = 2*pi*ones(block_height,1);
    end
end

%p(:,1) = 0.0005*ones(image_height,1);
%p(:,2) = -0.0314*ones(image_height,1);
%phase_0 = 63*ones(image_height,1);
%{        
for i = 1:image_height
    if p(i,1) > 0 
        %%%%%%%%%%%%%%%%%%位相接続のために絶対位相を二次関数で近似
        if phase_0(i) ~= 1
            for j = 1:phase_0(i)-1
                if theta(i,j) == 0
                    unwrapped_phase_temp(i,j) = 0;
                else 
                    unwrapped_phase_temp(i,j) = 0.5*p(i,1)*(phase_0(i)-j)^2+mean(theta(:,phase_0(i)));
                end
            end
        end
        if phase_0(i) ~= image_width
            for j = phase_0(i)+1:image_width
                if theta(i,j) ~= 0
                    unwrapped_phase_temp(i,j) = 0.5*p(i,1)*(j-phase_0(i))^2+mean(theta(:,phase_0(i)));
                else unwrapped_phase_temp(i,j) = 0;
                end
            end
        end
    else
        unwrapped_phase_temp(i,:) = zeros(1,image_width);
    end
end
%}

for i = 1:partition
    if p(i*block_height,1) > 0 
        %%%%%%%%%%%%%%%%%%位相接続のために絶対位相を二次関数で近似
        for k = 1:block_height
            for j = 1:image_width
                if theta(k+(i-1)*block_height,j) == 0
                    unwrapped_phase_temp(k+(i-1)*block_height,j) = 0;
                else 
                    unwrapped_phase_temp(k+(i-1)*block_height,j) = 0.5*p(i*block_height,1)*(phase_0(i*block_height)-j)^2+0.5262;
                end
            end
        end
    else
        unwrapped_phase_temp(i,:) = zeros(1,image_width);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%位相接続%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_map = medfilt2(unwrapped_phase_temp-theta,medfilt*2);
temp_out = temp_map/(2*pi);
index_temp = round(temp_map/(2*pi)); %位相接続の際のインデックス(仮)を格納
for i=1:image_height
    a = find(index_temp(i,:)-temp_out(i,:)~=0);
    if numel(a)~=0
        temp_out(i,:) = temp_out(i,:)+(index_temp(i,a(1))-temp_out(i,a(1)))*ones(1,image_width);
    end
end
index = round(temp_out); %位相接続の際のインデックスを格納

%{
for i = 1:image_height
    temp_map(i,:) = medfilt1(temp_map(i,:),10);
    temp_out(i,:) = temp_map(i,:)/gap(i);
    index(i,:) = round(temp_map(i,:)/gap(i));
end
%}




%{
for i= 1:image_height
    for j = 1:image_width-1
        if index(i,j)>index(i,j+1)
            index(i,j+1) = index(i,j);
        end
    end
end
%}

unwrapped_phase = index*2*pi + theta; %絶対位相を格納
unwrapped_phase = medfilt2(unwrapped_phase, medfilt);
          
%%%%%%%%%%%%%%%%%%%%%二次関数をプロット
%%グラフ5
figure;
plot(1:image_width,theta(line_out,:));
hold on;
plot(1:image_width,unwrapped_phase_temp(line_out,:),'r');
title('位相','FontSize',16);
xlim([0,image_width]);
set(gca,'FontSize',16);
saveas(gcf,['wrap_and_unwrap_',num2str(line_out),'.fig'])

%%%%%%%%%%%%%%%%%%%%%絶対位相をプロット
%%グラフ6
figure;
theta_lim = [-pi 5*2*pi];
unwrapped_phase = medfilt2(unwrapped_phase, medfilt);
imagesc(unwrapped_phase,theta_lim);
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
title('index','Fontsize',16);
saveas(gcf,['index_',num2str(line_out),'.fig'])


%%%%%%%%%%%%%%%%%%%%%%%%従来の位相接続
index2 = zeros(image_height,image_width);
for i = 1:image_height
    for j = 1:image_width-1
        if differ(i,j)< -pi
            index2(i,j+1:image_width) = index2(i,j+1:image_width) + ones(1,image_width-j);
        elseif theta(i,j)==0
            index2(i,j)=0;
        end
    end
end
unwrapped_phase2 = theta + index2*2*pi;
unwrapped_phase2 = medfilt2(unwrapped_phase2, medfilt);

%%%%従来の位相接続を補正するならこんな感じ
%{
for i = 1:image_height
    for j = 1:image_width
       if abs(unwrapped_phase2(i,j)-unwrapped_phase_temp(i,j))>2*pi
           unwrapped_phase2(i,j) = unwrapped_phase(i,j);
       end
    end
end
%}

%%%%%%従来の位相接続手法で求めた絶対位相をプロット
%%グラフ8
figure;
imagesc(unwrapped_phase2,theta_lim);
colorbar;
hold on;
plot(1:image_width,line_out,'r');
title('絶対位相（従来法）','Fontsize',16);
set(gca,'FontSize',16);
saveas(gcf,'unwrappedphase_2.fig')

cd('C:\Users\t2ladmin\Documents\MATLAB\phase');