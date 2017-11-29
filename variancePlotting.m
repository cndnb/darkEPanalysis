


  fixError = removeHighError(driftFix(:,1),driftFix(:,2),frequencyVariance(driftFix,1e-2,50),100);

  plot(fixError(:,1),fixError(:,3));
