[p_old,theta,differ,line_out,name] = plot_phase_old('\1125\hand_static - コピー',50);
[image_height,image_width] = size(theta);

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
%プロジェクタとの対応をとる
unwrapped_phase2(theta==-10) = NaN;
pixel2 = unwrapped_phase*64/pi;


%%%%%%従来の位相接続手法で求めた絶対位相をプロット
%%%グラフ
figure;
pixel2 = medfilt2(pixel2, medfilt);
imagesc(pixel2,pixel_lim);
colorbar;
hold on;
plot(1:image_width,line_out,'r');
title('絶対位相（従来法）','Fontsize',16);
set(gca,'FontSize',16);
saveas(gcf,'unwrappedphase_2.fig')
