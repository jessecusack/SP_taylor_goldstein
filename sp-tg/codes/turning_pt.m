%% This code pick out the turning points
% INPUTS:
% ww: profile, must be real or the imaginary part of a complex profile
%
% OUTPUTS:
% I: numbers of turning points
% l: indices
% 
% TO NOTICE:
%
% S.Tan, IOCAS, 2019/09/06
function [I,l]=turning_pt(ww,ref)

    I=nan(length(ww(1,:)),1);
    for k=1:length(ww(1,:))
        a=ww(:,k)-ref;
        aa=a(2:end).*a(1:end-1);
        l=find(aa<0);
        if isempty(l)
            l=0;
            I(k)=l;
        else
            I(k)=length(l);
        end
    end
