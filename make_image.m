function make_image()
    T=8;
    num = 1024/(2*T);
    A = zeros(768,1024);
    B = zeros(768,1024);
    C = zeros(768,1024);
    for t = 1: 1024
        %A(:,t)= (1+sin(pi/num*t-pi*2/3))/2;
        %B(:,t)= (1+sin(pi/num*t))/2;
        %C(:,t)= (1+sin(pi/num*t+pi*2/3))/2;
        A(:,t)= (1+sin(t^2/40000 - 2*pi/3))/2;
        B(:,t)= (1+sin(t^2/40000))/2;
        C(:,t)= (1+sin(t^2/40000 + 2*pi/3))/2;
        %A(:,t)= (1+sin(pi*(t-512)^2/40000 - 2*pi/3))/2;
        %B(:,t)= (1+sin(pi*(t-512)^2/40000 ))/2;
        %C(:,t)= (1+sin(pi*(t-512)^2/40000 + 2*pi/3))/2;
    end
    D = ones(768,1024)*0.5;
    str = 'sinimage_%d.bmp';
    imwrite(A,sprintf(str,1));
    imwrite(B,sprintf(str,2));
    imwrite(C,sprintf(str,3)); 
    imwrite(D,sprintf('reference.bmp')); 
end