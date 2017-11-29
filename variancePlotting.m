


  fixError = removeHighError(driftFix(:,1),driftFix(:,2),frequencyVariance(driftFix,1.333e-3,50),100);

  plot(fixError(:,1),fixError(:,3));
