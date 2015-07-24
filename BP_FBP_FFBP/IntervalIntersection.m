%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERVAL INTERSECTION FOR CALCULATING THE CONVOLUTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [IInt,RelAInd,RelBInd] = IntervalIntersection(IA,IB)

if IA(2) < IB(1) || IB(2) < IA(1) % --- No overlap case.
    IInt = [];
    RelAInd = [];
    RelBInd = [];
    return
end

if IB(2) < IA(2)
    IInt = [IA(1),IB(2)];
    RelAInd = 0;
    RelBInd = IA(1)-IB(1);    
    return
end

IInt = [IB(1),IA(2)];              % --- There is intersection to the right.
RelAInd = IB(1)-IA(1);
RelBInd = 0;
