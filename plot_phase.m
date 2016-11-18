%�擾�������摜������t�H���_������͂���
%�ʑ��̉�A�W���C�e�s�N�Z���̈ʑ��Cx�����̈ʑ��̔������o��
function [p,theta,differ,line_out,name] = plot_phase(name,line_out)
    close all;
    %%%%%%%%%%
    %%%%%%%%%%
    
    folder_name = 'C:\Users\t2ladmin\Documents\MATLAB\phase shift\phase_image';
    file_name = [folder_name, '\', name];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�t�H���_���Ɉړ����C�摜�݂̂�ǂݍ���
    cd(file_name);
    imgFileList = dir('*.tif'); 
    
    num = length(imgFileList); %�摜�̐�
    image_width = 512; %�J�����̃s�N�Z����(x����)
    image_height = 512; %�J�����̃s�N�Z����(y����)
    pixel = image_width * image_height; %�J�����̑S�s�N�Z����
    Luminance = zeros(pixel,num); %�S�摜�C�S�s�N�Z���̋P�x���i�[����s��

    %���ۂɋP�x���L�^����
    for i = 1:num
        img_name = char(imgFileList(i).name);
        temp = imread(img_name);  %�J�����̎B�e�摜�̓ǂݍ���
        %�J���[�摜�̏ꍇ2�l������
        if(numel(temp) == pixel*3)
            temp = rgb2gray(temp);
        end
        temp = temp(:);
        Luminance(:,i) = temp;
        clear temp
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�P�x�̋L�^�I��
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�ʑ��̌v�Z
    Luminance = Luminance';
    M = [1,-1/2,sqrt(3)/2;1,1,0;1,-1/2,-sqrt(3)/2];
    temp = M\Luminance;
    clear M;
    I_mod = sqrt(temp(2,:).^2+temp(3,:).^2); %�U���̓��o
   
    theta = zeros(1,pixel); %�ʑ����i�[
    %�ʑ����v�Z���Ă���
    for i = 1:pixel
        if I_mod(i)<5
            theta(i)=0; %�U�����������Ƃ��C���m�Ȉʑ��̌v�Z�͕s�\
        else
            cos = temp(2,i)./I_mod(i);
            sin = temp(3,i)./I_mod(i);
            if abs(cos)<=1 && abs(sin)<=1
                theta(i) = atan2(sin,cos);
            else
                theta(i)=0;
            end
        end
    end
    
    theta = reshape(theta,[image_height image_width]); %-pi����pi�͈̔͂ňʑ������߂�
    medfilt = [5 5];
    theta = medfilt2(theta, medfilt);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�ʑ��̌v�Z�I��
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�ʑ��̔����ɒ���
    differ = zeros(image_height,image_width-1); %theta�̍��ix�����̔����j���i�[
    p = zeros(image_height,2); %theta�̔����̉�A���s���ƂɊi�[
    for line = 1:image_height
        %�ʑ��̋�ԕ����̔������v�Z
        differ(line,:) = diff(theta(line,:));
    
        %�O��l�����O���C�ʑ��̔�������`�ߎ�
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�ʑ��̔����̌v�Z�Ɖ�A�I��
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����s�ɂ��Ă̈ʑ���P�x���v���b�g����
    %%line_out = 150;
    %�ʑ��̃v���b�g
    figure;
    plot(1:image_width,theta(line_out,:));
    title('�ʑ�','FontSize',16);
    xlim([0,image_width]);
    ylim([-pi,pi]);
    set(gca,'FontSize',16);
    saveas(gcf,'phase.fig')
    
    %�ʑ��̋�ԕ����̔������v���b�g
    figure;
    plot(1:image_width-1,differ(line_out,:));
    hold on
    plot(1:image_width-1,polyval(p(line_out,:),1:image_width-1),'r'); %�ʑ��̔����̐��`�ߎ����v���b�g
    title('�ʑ��̔���','FontSize',16);
    xlim([0,image_width]);
    ylim([-0.5,0.5]);
    set(gca,'FontSize',16);
    saveas(gcf,'derivative.fig')
    
    %�P�x���v���b�g�Csin�g�ɂȂ��Ă��OK
    figure;
    sin_wave = reshape(Luminance(2,:), [image_height image_width]);
    plot(1:image_width,sin_wave(line_out,:));
    title('�P�x','FontSize',16);
    xlim([0,image_width]);
    ylim([0,256]);
    set(gca,'FontSize',16);
    saveas(gcf,'I.fig')
  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ۂ̈ʑ������摜��\���D�s�v�ȏꍇ�͈ȉ����R�����g�A�E�g
    figure;
    %�����Ȃ�قǈʑ���2pi�ɋ߂Â�
    theta_lim = [-pi pi];
    imagesc(theta,theta_lim);
    colorbar;
    hold on;
    plot(1:image_width,line_out,'r');
    title('�ʑ��}�b�v','FontSize',16);
    set(gca,'FontSize',16);
    saveas(gcf,'phasemap.fig')
    
    %%���̃f�B���N�g���ɖ߂�
    cd('C:\Users\t2ladmin\Documents\MATLAB\phase shift');
end