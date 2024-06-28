function [rmserr,rmsalt] = rmsecalc(fitval,dataval,numcoeff)
% RMS deviation
%     numval = 0;
%     sumerr = 0;
%     for i = 1:length(dataval)
%         if ~isnan(dataval(i))
%             sumerr = sumerr + (fitval(i)-dataval(i))^2;
%             numval = numval + 1;
%         end
%     end
%     if (numval > 0)
%         rmserr = sqrt(sumerr/numval)
%     else
%         rmserr = 0;
%     end
%     if isnan(sumerr)
%         [fitval dataval]'
%     end
    rmsindex = ~isnan(dataval);         % Indices of valid data values
    fitvalrms = fitval(rmsindex);       % Fit at existing data points
    datavalrms = dataval(rmsindex);     % Data at existing data points
    % val = [fitval dataval]';
    % valrms = [fitvalrms datavalrms]';
    numval = sum(rmsindex);             % Number of valid data values
    % sumerr = sum((fitval - dataval).^2) does not account for missing data
    sumerr = sum((fitvalrms - datavalrms).^2);
    rmserr = sqrt(sumerr/numval);
    dfe = numval - numcoeff;
    if (dfe > 0)
        rmsalt = sqrt(sumerr/dfe);
    else
        rmsalt = NaN;
    end
end
