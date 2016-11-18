[p,theta,differ,line_out,name] = plot_phase('\1103_experiment_contrast\plane_2000_cont20_mono',350);
[image_height,image_width] = size(theta);
unwrapped_phase_temp = zeros(image_height,image_width);
index = zeros(image_height,image_width);
unwrapped_phase = zeros(image_height,image_width);

folder_name = 'C:\Users\t2ladmin\Documents\MATLAB\phase shift\phase_image';
file_name = [folder_name, '\', name];
cd(file_name);

%%絶対位相が0となるピクセルを推定
temp_p1= p(:,1);
A = find(temp_p1<=0);
temp_p2= p(:,2);
B = find(temp_p2>=0);
temp_p = p;
for i = 1:length(A)
    temp_p(A(i),:)=[0,0];
end
for i = 1:length(B)
    temp_p(B(i),:)=[0,0];
end
temp_p(temp_p==0)=[];
temp_p = reshape(temp_p,[length(temp_p)/2 2]);
phase_0 = round(-mean(temp_p(:,2))/mean(temp_p(:,1)))*ones(image_height,1);
p = horzcat(mean(temp_p(:,1))*ones(image_height,1),mean(temp_p(:,2))*ones(image_height,1));

%{
phase_0 = zeros(image_height,1);
for i = 1:image_height
    if p(i,1) > 0
        %%絶対位相が0となるピクセルを推定
        phase_0(i) = round(-p(i,2)/p(i,1));
        if phase_0(i)<1
            phase_0(i) = 1;
        elseif phase_0(i)>512
            phase_0(i) = 1;
        end
    end
end
plot(phase_0,1:image_width);
%}

        
for i = 1:image_height
    if p(i,1) > 0 
        %%位相接続
        if phase_0(i) ~= 1
            for j = 1:phase_0(i)-1
                if theta(i,j) == 0
                    unwrapped_phase_temp(i,j) = 0;
                else 
                    unwrapped_phase_temp(i,j) = -0.5*p(i,1)*(phase_0(i)-j)^2+mean(theta(:,1));
                end
            end
        end
        if phase_0(i) ~= image_width
            for j = phase_0(i)+1:image_width
                if theta(i,j) ~= 0
                    unwrapped_phase_temp(i,j) = 0.5*p(i,1)*(j-phase_0(i))^2+mean(theta(:,1));
                else unwrapped_phase_temp(i,j) = 0;
                end
            end
        end
    else
        unwrapped_phase_temp(i,:) = zeros(1,image_width);
    end
    
    index = round((unwrapped_phase_temp-theta)/(2*pi));
    for i = 1:image_height
        for j = 1:image_width
            if j<phase_0(i)
                index(i,j) = -index(i,j);
                theta(i,j) = -theta(i,j);
            end
        end
    end
    unwrapped_phase = index*2*pi + theta;
end
          
%{        
for i = 1:image_height
    if p(i,1) > 0
        for j = 1:image_width
            if theta(i,j) > -3.5
                unwrapped_phase(i,j) = theta(i,j) + 0.5*p(i,1)+j^2 + p(i,2)*j;
            else
                unwrapped_phase(i,j) = -3.5;
            end
        end
    end
end
%}

%%%%%%%%絶対位相をプロット
figure;
theta_lim = [-pi 7.5*2*pi];
medfilt = [5 5];
unwrapped_phase = medfilt2(unwrapped_phase, medfilt);
imagesc(unwrapped_phase,theta_lim);
colorbar;
%{
theta_img = mat2gray(unwrapped_phase);
imshow(theta_img);
%}
hold on;
plot(1:image_width,line_out,'r');
title('絶対位相','Fontsize',16);
set(gca,'FontSize',16);
saveas(gcf,'unwrappedphase_1.fig')

%%%%%%位相接続時のインデックスをプロット
figure;
plot(1:image_width,index(line_out,:));
xlim([0,image_width]);
title('index','Fontsize',16);



%%%%%%従来の位相接続
index2 = zeros(image_height,image_width);
for i = 1:image_height
    for j = 1:image_width-1
        if differ(i,j)< -pi
            index2(i,j+1:image_width) = index2(i,j+1:image_width) + ones(1,image_width-j);
        end
    end
end
unwrapped_phase2 = theta + index2*2*pi;
%%%%%%従来の位相接続手法で求めた絶対位相をプロット
figure;
imagesc(unwrapped_phase2,theta_lim);
colorbar;
hold on;
plot(1:image_width,line_out,'r');
title('絶対位相（従来法）','Fontsize',16);
set(gca,'FontSize',16);
saveas(gcf,'unwrappedphase_2.fig')

cd('C:\Users\t2ladmin\Documents\MATLAB\phase shift');