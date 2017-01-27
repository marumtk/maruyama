[p_old,theta,differ,line_out,name] = plot_phase_old('\1125\hato_static',200);
[image_height,image_width] = size(theta);
unwrapped_phase_temp = zeros(image_height,image_width);
medfilt = [5 5];


folder_name = 'C:\Users\t2ladmin\Documents\MATLAB\phase shift\phase_image';
file_name = [folder_name, '\', name];
cd(file_name);

%%%%%%%%%%%%%%%%%%%%%%���������i�[,�O��l������
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
%%%%%%%%%%%%%%%%%%%%%%%%��Έʑ���0�ƂȂ�s�N�Z���𐄒�
phase_0 = round(-mean(temp_p(:,2))/mean(temp_p(:,1)))*ones(image_height,1);
p = horzcat(mean(temp_p(:,1))*ones(image_height,1),mean(temp_p(:,2))*ones(image_height,1));
%phase_0 = 132*ones(image_height,1);
%p = horzcat(0.0007*ones(image_height,1),-0.0926*ones(image_height,1));

%%%%%%%%%%%%%%%%%�ʑ��̔��������v���b�g
%�O���t4
figure;
plot(1:image_width-1,differ(line_out,:));
hold on
plot(1:image_width-1,polyval(p(line_out,:),1:image_width-1),'r'); %�ʑ��̔����̐��`�ߎ����v���b�g
title(['�ʑ��̔���',num2str(line_out)],'FontSize',16);
xlim([0,image_width]);
ylim([-0.5,0.5]);
set(gca,'FontSize',16);
saveas(gcf,['derivative_',num2str(line_out),'.fig'])


%%%%%%%%%%%%%%%%%������񂩂�ʑ���\��2���֐��𓱏o       
for i = 1:image_height
    if p(i,1) > 0 
        %%�ʑ��ڑ�
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
          
%%%%%%%%%%%%%%%%%%%%%%%%%�ʑ��ڑ�%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_map = unwrapped_phase_temp-theta;
temp_out = temp_map/(2*pi);
index = round(temp_out); %�ʑ��ڑ��̍ۂ̃C���f�b�N�X���i�[
unwrapped_phase = index*2*pi + theta; %��Έʑ����i�[
%�v���W�F�N�^�Ƃ̑Ή����Ƃ�
pixel = 200*sqrt(abs(unwrapped_phase));

%%%%%%%%%%%%%%%%%%%%%�񎟊֐����v���b�g
%%�O���t5
figure;
plot(1:image_width,theta(line_out,:));
hold on;
plot(1:image_width,unwrapped_phase_temp(line_out,:),'r');
title(['�ʑ�',num2str(line_out)],'FontSize',16);
xlim([0,image_width]);
set(gca,'FontSize',16);
saveas(gcf,['wrap_and_unwrap_',num2str(line_out),'.fig'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��Έʑ����v���b�g
%%�O���t6
figure;
pixel_lim = [0 1024];
%pixel = medfilt2(pixel, medfilt);
imagesc(pixel,pixel_lim);
colorbar;
hold on;
plot(1:image_width,line_out,'r');
title('��Έʑ�','Fontsize',16);
set(gca,'FontSize',16);
saveas(gcf,'unwrappedphase_1.fig')

%%%%%%%%%%%%%%%%%%%%%�ʑ��ڑ����̊(�񎟊֐�)�ƃ��b�v���ꂽ�ʑ��̍����v���b�g
%%�O���t7
figure;
plot(1:image_width,temp_out(line_out,:),'r');
hold on;
plot(1:image_width,index(line_out,:),'g');
hold on;
plot(1:image_width,theta(line_out,:));
xlim([0,image_width]);
title(['index',num2str(line_out)],'Fontsize',16);
saveas(gcf,['index_',num2str(line_out),'.fig'])

%%%%%%�]���̈ʑ��ڑ�
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
%%%%%%�]���̈ʑ��ڑ���@�ŋ��߂���Έʑ����v���b�g
%%%�O���t8
figure;
imagesc(pixel2,pixel_lim);
colorbar;
hold on;
plot(1:image_width,line_out,'r');
title('��Έʑ��i�]���@�j','Fontsize',16);
set(gca,'FontSize',16);
saveas(gcf,'unwrappedphase_2.fig')

figure
plot(1:image_width,pixel(line_out,:));
title('�v���W�F�N�^�s�N�Z��','Fontsize',16);
xlabel('�J������','Fontsize',16)
ylabel('�s�N�Z��','Fontsize',16)
ylim([0 1024]);
set(gca,'FontSize',16);
saveas(gcf,['pixel_info',num2str(line_out),'.fig'])



cd('C:\Users\t2ladmin\Documents\MATLAB\phase shift');