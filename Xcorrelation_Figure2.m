function XCorr = Xcorrelation_Figure5(cTrials, LagSteps, SRAnalysis)


NTrials = numel(cTrials);


 for iT1 = 1:NTrials-1
      cTimes1 = cTrials{iT1};
      cIndices1 = round(cTimes1/SRAnalysis);
      for iT2 = iT1+1:NTrials
        cTimes2 = cTrials{iT2};
        cIndices2 = round(cTimes2/SRAnalysis);
        cX = xcorrtimes(cIndices1,cIndices2,LagSteps,1);
        XCorr(:) = XCorr(:) + cX';
      end
 end



end




function out = xcorrtimes(in1,in2,maxlag,factor)


if (nargin < 4) factor = 1; end
if (nargin < 3) maxlag = 100; end
if (nargin < 2) in2 = in1; end

if factor ~=1 
  in1 = res2res(in1,factor);
  in2 = res2res(in2,factor);
end

maxlength = max([in1,in2]);
in1_time = time2event(in1,maxlength);  
in2_time = time2event(in2,maxlength);

out = binxcorr(in1_time',in2_time',maxlag,1);
end

function out = res2res(in,factor)
% out = res2res(in,factor) : change resolution
% in : vector containing values in first resolution
% factor : change in resolution
% out : vector containing values in second resolution

out = in*factor;
end

function out = time2event(in,tl)
% time2event(in,tl)
% in : vector containing event-times
% tl : specifies the total length of the vector
% converts eventtimes into a vector, where the position in the vector 
% signifies the time of the event.

out = zeros(tl,1);
in = round(in);
i=1;
while (i <=length(in) & in(i) <= tl)
  if (in(i) ~= 0) out(in(i),1) = out(in(i),1) + 1; end;
  i = i + 1;
end

end


function out = binXcorr(in1,in2,maxlag,option)
% binXcorr(in1,in2,maxlag,option)
% in1 : vector of events
% in2 : vector of events or cell of times : depends on option
% maxlag : maximal lag between the two sequences analyzed
% option : 1: in2 is a matrix of events, 2: in2 is a cell of times
%
% binXcorr is a fast replacement for xcorr, if only sparse, binary sequences 
% are to be analyzed. 
% 
% see also: ACXCORR, WINXCORR, XCORR

if (nargin < 4) option = 1; end    
if (nargin < 3) maxlag = 100; end
if (nargin < 2) in2 = in1; end
 
out = zeros(1,2*maxlag+1);   

if option == 1   % matrix and matrix
    in1 = [zeros(1,maxlag),in1,zeros(1,maxlag)];
    in2 = [zeros(1,maxlag),in2,zeros(1,maxlag)];
    if (length(in1)<length(in2)) arg = length(in1); else arg = length(in2); end; 
    for i=maxlag+1:arg-maxlag-1
        if in1(i) == 1
            out = out + in2(i-maxlag:i+maxlag);
        end
    end
    
elseif option == 2   % cell and matrix 
    in2 = in2 + maxlag;
    in1 = double([zeros(1,maxlag),in1,zeros(1,maxlag)]);
    for i = 1:length(in2)
        out = out + in1(in2(i)-maxlag:in2(i)+maxlag,:);
    end    
    
else error('option not recognized: 1 for binary, 2 for times')
end
out = fliplr(out);
end