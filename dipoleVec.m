%Constructs a time dependent Seattle dipole vector
%Arguments in units of (seconds, decimal degrees, " ", " ", magnitude of dipole)
function ret = dipoleVec(t, seattleLat, seattleLong, compassDir, dM)
  global omegaEarth;
  A = ones(length(t),1).*dM*cos(deg2rad(compassDir))*cos(deg2rad(seattleLong));
  B = dM*sin(deg2rad(compassDir))*sin(deg2rad(seattleLat) + omegaEarth.*t);
  C = dM*sin(deg2rad(compassDir))*cos(deg2rad(seattleLat) + omegaEarth.*t);
  ret = [t, A, B, C];
endfunction