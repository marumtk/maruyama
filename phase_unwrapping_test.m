[p_old,theta,differ,line_out,name] = plot_phase_old('\1220\hand_proposed',350);
[image_height,image_width] = size(theta);
unwrapped_phase_temp = zeros(image_height,image_width);
phase_0 = zeros(image_height,1); %�ʑ��̔�����0�ƂȂ�_�𐄒�
p = zeros(image_height,2); %���ϔ��������i�[
medfilt = [10 10];


folder_name = 'C:\Users\t2ladmin\Documents\MATLAB\phase shift\phase_image';
file_name = [folder_name, '\', name];
cd(file_name);

%%%%%%%%%%%%%%%%%%%%%%���������i�[,�O��l������
temp_p1= p_old(:,1);
A = find(temp_p1<=0);
temp_p = p_old;
Ind = 1:768; %%��������1024*768�̉摜��p���鎞��Ind�����g���ς��

%%���̉摜�f�[�^��p���鎞�͈ȉ����R�����g�A�E�g
%�O��l����
for i = 1:length(A)
    temp_p(A(i),:)=[NaN,NaN];
end

%%�s�N�Z���������Ȃ��s������
for i = 1:512
    Temp = find(theta(i,:)~=-10);
    if length(Temp) < 50
        temp_p(i,:) = [NaN,NaN];
    end
end
Ind = find(isnan(temp_p(:,1))==0); %�l�������Ă���s���i�[
temp_p(isnan(temp_p)==1)=[];
temp_p = reshape(temp_p,[length(temp_p)/2 2]);

%�[������
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
%%���̉摜�f�[�^��p���鎞�͂����܂ł��R�����g�A�E�g

%%%%%%%%%%%%%%%%%%%%%%%%��Έʑ���0�ƂȂ�s�N�Z���𐄒�
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
          
%%%%%%%%%%%%%%%%%%%%%%%%%�ʑ��ڑ�%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_map = unwrapped_phase_temp-theta;
temp_out = temp_map/(2*pi);
index = round(temp_out); %�ʑ��ڑ��̍ۂ̃C���f�b�N�X���i�[
unwrapped_phase = index*2*pi + theta; %��Έʑ����i�[
unwrapped_phase = medfilt2(unwrapped_phase,medfilt);
unwrapped_phase(unwrapped_phase==-10) = NaN;
%�v���W�F�N�^�Ƃ̑Ή����Ƃ�
pixel = 200*sqrt(abs(unwrapped_phase));
%�덷���傫���Ƃ��͓񎟊֐��𗘗p
%pixel_temp = 200*sqrt(abs(unwrapped_phase_temp));
%pixel(abs(pixel_temp-pixel)>1) = pixel_temp(abs(pixel_temp-pixel)>1);

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
colormap(jet(256));
axis image;
hold on;
%plot(1:image_width,line_out,'r');
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

%%%%%%%%%%%%%%%%%%%%%�v���W�F�N�^�̃s�N�Z�����v���b�g
%%�O���t8
figure
plot(1:image_width,pixel(line_out,:));
title('�v���W�F�N�^�s�N�Z��','Fontsize',16);
xlabel('�J������','Fontsize',16)
ylabel('�s�N�Z��','Fontsize',16)
xlim([0 image_width]);
set(gca,'FontSize',16);
saveas(gcf,['pixel_info',num2str(line_out),'.fig'])


cd('C:\Users\t2ladmin\Documents\MATLAB\phase shift');