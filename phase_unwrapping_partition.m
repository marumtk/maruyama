[p,theta,differ,line_out,name,gap] = plot_phase('\1125\hand_static_noref',50);
[image_height,image_width] = size(theta); %�J�����s�N�Z���̃T�C�Y

unwrapped_phase_temp = zeros(image_height,image_width); %���̐�Έʑ����i�[(�񎟊֐��ߎ��̌���)
partition = 1; %�����������߂�ۂ̕�����(�u���b�N��)
block_height = image_height/partition; %�����������߂�ۂ̈�u���b�N�̃T�C�Y
partition_noise = zeros(partition,2); %��u���b�N������̃m�C�Y��
medfilt = [5 5]; %���f�B�A���t�B���^

folder_name = 'C:\Users\t2ladmin\Documents\MATLAB\phase\phase_image';
file_name = [folder_name, '\', name];
cd(file_name);


%%%%%%%%%%%%%%%%%%%%%%���������i�[,�O��l������
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


%%%%%%%%%%%%%%%%%%%%%�񎟊֐��̒��_�ƂȂ�s�N�Z����񂲂Ƃɐ���
phase_0 = zeros(image_height,1); %theta^2��\���񎟊֐��̒��_�ƂȂ�s�N�Z�����i�[

%%%%%%%%%%%%����ʂ�p�ӂ��ĕ��s�X�e���I�ƂȂ�悤�ɎB������Ȃ炱�̃A���S���Y�����g����
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

%%%%%%%%%%%%%%%%%�m�C�Y�̐����u���b�N���ƂɎZ�o
for i = 1:partition
    num1 = find(temp_p(1+(i-1)*block_height:i*block_height,1)==0);
    num2 = find(gap(1+(i-1)*block_height:i*block_height,1)==0);
    partition_noise(i,1) = length(num1); %�������̃m�C�Y�̐�
    partition_noise(i,2) = length(num2); %�񎟊֐��̒��_�ɂ��Ẵm�C�Y�̐�
end

%%%%%%%%%%%%%%%%%�������̕��ϒl���u���b�N���Ƃɓ��o        
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
        %%%%%%%%%%%%%%%%%%�ʑ��ڑ��̂��߂ɐ�Έʑ���񎟊֐��ŋߎ�
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
        %%%%%%%%%%%%%%%%%%�ʑ��ڑ��̂��߂ɐ�Έʑ���񎟊֐��ŋߎ�
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

%%%%%%%%%%%%%%%%%%%%%%%%%�ʑ��ڑ�%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_map = medfilt2(unwrapped_phase_temp-theta,medfilt*2);
temp_out = temp_map/(2*pi);
index_temp = round(temp_map/(2*pi)); %�ʑ��ڑ��̍ۂ̃C���f�b�N�X(��)���i�[
for i=1:image_height
    a = find(index_temp(i,:)-temp_out(i,:)~=0);
    if numel(a)~=0
        temp_out(i,:) = temp_out(i,:)+(index_temp(i,a(1))-temp_out(i,a(1)))*ones(1,image_width);
    end
end
index = round(temp_out); %�ʑ��ڑ��̍ۂ̃C���f�b�N�X���i�[

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

unwrapped_phase = index*2*pi + theta; %��Έʑ����i�[
unwrapped_phase = medfilt2(unwrapped_phase, medfilt);
          
%%%%%%%%%%%%%%%%%%%%%�񎟊֐����v���b�g
%%�O���t5
figure;
plot(1:image_width,theta(line_out,:));
hold on;
plot(1:image_width,unwrapped_phase_temp(line_out,:),'r');
title('�ʑ�','FontSize',16);
xlim([0,image_width]);
set(gca,'FontSize',16);
saveas(gcf,['wrap_and_unwrap_',num2str(line_out),'.fig'])

%%%%%%%%%%%%%%%%%%%%%��Έʑ����v���b�g
%%�O���t6
figure;
theta_lim = [-pi 5*2*pi];
unwrapped_phase = medfilt2(unwrapped_phase, medfilt);
imagesc(unwrapped_phase,theta_lim);
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
title('index','Fontsize',16);
saveas(gcf,['index_',num2str(line_out),'.fig'])


%%%%%%%%%%%%%%%%%%%%%%%%�]���̈ʑ��ڑ�
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

%%%%�]���̈ʑ��ڑ���␳����Ȃ炱��Ȋ���
%{
for i = 1:image_height
    for j = 1:image_width
       if abs(unwrapped_phase2(i,j)-unwrapped_phase_temp(i,j))>2*pi
           unwrapped_phase2(i,j) = unwrapped_phase(i,j);
       end
    end
end
%}

%%%%%%�]���̈ʑ��ڑ���@�ŋ��߂���Έʑ����v���b�g
%%�O���t8
figure;
imagesc(unwrapped_phase2,theta_lim);
colorbar;
hold on;
plot(1:image_width,line_out,'r');
title('��Έʑ��i�]���@�j','Fontsize',16);
set(gca,'FontSize',16);
saveas(gcf,'unwrappedphase_2.fig')

cd('C:\Users\t2ladmin\Documents\MATLAB\phase');