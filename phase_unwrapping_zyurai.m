[p_old,theta,differ,line_out,name] = plot_phase_old('\1125\hand_static - �R�s�[',50);
[image_height,image_width] = size(theta);

%%%%%%�]���̈ʑ��ڑ�
index2 = zeros(image_height,image_width);
for i = 1:image_height
    for j = 1:image_width-1
        if differ(i,j)< -pi
            index2(i,j+1:image_width) = index2(i,j+1:image_width) + ones(1,image_width-j);
        end
    end
end
unwrapped_phase2 = theta + index2*2*pi;
%�v���W�F�N�^�Ƃ̑Ή����Ƃ�
unwrapped_phase2(theta==-10) = NaN;
pixel2 = unwrapped_phase*64/pi;


%%%%%%�]���̈ʑ��ڑ���@�ŋ��߂���Έʑ����v���b�g
%%%�O���t
figure;
pixel2 = medfilt2(pixel2, medfilt);
imagesc(pixel2,pixel_lim);
colorbar;
hold on;
plot(1:image_width,line_out,'r');
title('��Έʑ��i�]���@�j','Fontsize',16);
set(gca,'FontSize',16);
saveas(gcf,'unwrappedphase_2.fig')
