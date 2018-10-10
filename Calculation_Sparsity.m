function S = calcSparsity_Figure8(cPSTH,Method)


switch Method

    case 'CloseToSilent' 
        Threshold = 15;
        S = sum(cPSTH<Threshold)/length(cPSTH);
    
    case 'Variance' 
        S = 1 - (sum(cPSTH)/length(cPSTH)).^2/(sum(cPSTH.^2)/length(cPSTH));
     
    case 'Kurtosis' 
        S = kurtosis(cPSTH);
        
      case 'DistanceToConst' % 
        cPSTH = cPSTH/max(cPSTH);
        S = sum(ones(size(cPSTH)) - cPSTH)/length(cPSTH);
   
    case 'SmallerThanMean'
        S = sum(cPSTH < mean(cPSTH))/length(cPSTH);
    
    case 'Skewness' 
        cPSTHN = (cPSTH - mean(cPSTH))/std(cPSTH);
        S = skewness(cPSTHN);
    
    otherwise error('Method not implemented.');
end


end