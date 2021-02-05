%butterworth low-pass filter ;xn is series
function [fff_low,H,w]=butterworth_lp(xn,fl,order,point)
dt=1;
fff=xn;
wn_low=1/fl/(1/(2*dt));%%%%%%f1=time,xn=series
%wn cutoff frequency
[filter_lowa,filter_lowb]=butter(order,wn_low,'low');
[H,w]=freqz(filter_lowa,filter_lowb,point);                   
fff_low=filtfilt(filter_lowa,filter_lowb,fff);
return



