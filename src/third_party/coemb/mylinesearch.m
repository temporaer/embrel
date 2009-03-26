function [xc,stepsize] =  mylinesearch(x0,f0,f,par,dsd,maxstep,stepfactor)

stepsize = maxstep;
fc = f0+1;

while fc > f0
  xc = x0 + stepsize*dsd;
  fc = feval(f,xc,par);
  stepsize = stepsize*stepfactor;
end

fprintf('Final factor %g\n',stepsize/stepfactor);
