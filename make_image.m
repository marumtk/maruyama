function make_image(T)
    num = 1024/(2*T);
    A = zeros(768,1024);
    B = zeros(768,1024);
    C = zeros(768,1024);
    for t = 1: 1024
        %A(:,t)= (1+sin(pi/num*t))/2;
        %B(:,t)= (1+sin(pi/num*t+pi*2/3))/2;
        %C(:,t)= (1+sin(pi/num*t+pi*4/3))/2;
        A(:,t)= (1+sin(pi*t^2/65536))/2;
        B(:,t)= (1+sin(pi*t^2/65536 + 2*pi/3))/2;
        C(:,t)= (1+sin(pi*t^2/65536 + 4*pi/3))/2;
        %A(:,t)= (1+cos(pi*(t-512)^2/40000))/2;
        %B(:,t)= (1+cos(pi*(t-512)^2/40000 + 2*pi/3))/2;
        %C(:,t)= (1+cos(pi*(t-512)^2/40000 + 4*pi/3))/2;
    end
    str = 'sinimage_%d.tif';
    imwrite(A,sprintf(str,1));
    imwrite(B,sprintf(str,2));
    imwrite(C,sprintf(str,3)); 
    imshow((A+B)/2);
end