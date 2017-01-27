function [out_low] = lowpass_filt(input,Fs)
    b_lp = fir1(1,Fs);
    out_low = filter(b_lp,1,input);
    t = (0:length(input)-1);
    plot(t,out_low);
    hold on;
    plot(t,input,'r');
end
