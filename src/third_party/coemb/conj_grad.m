function  [X, info, perf] = conj_grad(fun,par, x0, opts)
%CONJ_GRAD  Conjugate gradient method for nonlinear optimization:
%  Find  xm = argmin{f(x)} , where  x  is an n-vector and the scalar
%  function  f  with gradient  g  (with elements  g(i) = Df/Dx_i )
%  must be given by a MATLAB function with with declaration
%            function  [f, g] = fun(x, par)
%  par  holds parameters of the function.  It may be dummy.  
%  
%  Call
%      [X, info {, perf}] = conj_grad(fun,par, x0, opts)
%
%  Input parameters
%  fun  :  String with the name of the function.
%  par  :  Parameters of the function.  May be empty.
%  x0   :  Starting guess for  x .
%  opts :  Vector with 9 elements:
%          opts(1) :  Choice of CG method:
%              opts(1) = 1 :  Fletcher - Reeves, 
%                otherwise :  Polak - Ribiere 
%          opts(2) :  Choice of line search method:
%              opts(2) = 1 :  Exact l. s., otherwise : Soft l.s.
%          opts(3) :  Upper bound on initial step
%          opts(4:6)  used in stopping criteria:
%              ||g||_inf <= opts(4)                     or 
%              ||dx||_2 <= opts(5)*(opts(5) + ||x||_2)  or
%              no. of function evaluations exceeds  opts(6) . 
%          opts(7:9)  used as linesearch parameters lpar(3:5) .
%              See LINESEARCH .
%          Any illegal element in  opts  is replaced by its
%          default value :
%              opts(1:6) = [2  2  1  1e-4*||g(x0)||_inf  1e-6  100]
%          and
%              opts(2)=1:  opts(7:9) = [1e-4  1e-6  10]
%              opts(2)=2:  opts(7:9) = [1e-1  1e-2  10]
%
%  Output parameters
%  X    :  If  perf  is present, then array, holding the iterates
%          columnwise.  Otherwise, computed solution vector.
%  info :  Performance information, vector with 6 elements:
%          info(1:3) = final values of 
%              [f(x)  ||g||_inf  ||dx||_2] 
%          info(4:5) = no. of iteration steps and evaluations of (f,g)
%          info(6) = 1 :  Stopped by small gradient
%                    2 :  Stopped by small x-step
%                    3 :  Stopped by  opts(6) .  
%  perf :  (optional). If present, then array, holding 
%            perf(1:2,:) = values of  f(x) and ||g||_inf
%            perf(3:5,:) = Line search info:  values of  
%                          alpha, phi'(alpha), no. fct. evals.
%
%  Use of other MATLAB functions :  LINESEARCH .
%
%  Method
%  See Chapter 4 in  P.E. Frandsen, K. Jonasson, H.B. Nielsen, 
%  O. Tingleff: "Unconstrained Optimization", IMM, DTU.  1999.

   verbose = 1;
%  Hans Bruun Nielsen,  IMM, DTU.  99.08.05-09 / 08.23 / 12.07
    jitter_mag = 0.01;
    alfa=-1;
    if isfield(par,'verbose')
      verbose = par.verbose;
    end
   %  Check call 
   if length(opts)==10
     nMaxIter = opts(10);
     verbose = 1;
   else
     nMaxIter = 0;
   end
   opts = opts(1:9);
   [x  n  f  g  opts] = check(fun,par,x0,opts);

   %  Finish initialization
   k = 1;   kmax = opts(6);
   ng = norm(g, inf);   n2g = norm(g);   gg = n2g^2;
   gam = 0;   h = zeros(size(x));   nh = 0;   neval = 1;
   Trace = nargout > 2;
   if  Trace
         X = x(:)*ones(1,kmax+1);
         perf = [f; ng; zeros(3,1)]*ones(1,kmax+1);
       end
%   opts(4)=1e-12; %Amir
   found = ng <= opts(4);
   nit=0;
   
   while  ~found
     if verbose
         tic
     end
     %  Previous values
     hpr = h;   nhpr = nh;  ggpr = gg; nit=nit+1;
     if  opts(1) == 2,   gpr = g; end
     %  New direction and scaling
     k = k+1;   h = gam*h - g;   nh = norm(h);
     %  Check descent
     gh = dot(g,h);     
     if  dot(g,h) >= -1e-3 * n2g * nh    % Not descent
%       disp('Flipping Descent');
       h = -g;   nh = n2g;
     end
     % Scale, aiming at  alpha = 1
     if  k > 2,  sc = .9*alfa*nhpr/nh;
     else    % First step.  Try not to move too far
       sc = opts(3)/32/nh;
     end
     %  Line search
     
     lspar = [32  opts([2 7:8])  min(kmax-neval, opts(9))];

     [al2  f2  g2  dval2] = linesearch(fun,par,x,f,g, sc*h, lspar);
     if al2==0
       %print_to_log('Using it');
%       x = x + jitter_mag/nit*2*(0.5-rand(size(x))).*x;
%       lspar = [32  opts([2 7:8]) min(kmax-neval,opts(9))];
       h = -g;
       [al  f  g  dval] = linesearch(fun,par,x,f,g, sc*h, lspar); 
       %print_to_log('Now Alpha is %g\n',al);
     else
       al = al2; f = f2; g = g2; dval = dval2;
     end
	   
     alfa = sc*al;   ng = norm(g, inf);   n2g = norm(g);  gg = n2g^2;
     %  Update  gamma  and  x
     if  opts(1) == 1,  gam = gg/ggpr;
     else,              gam = dot(g-gpr,g)/ggpr; end
     x = x + alfa*h;  
%     x = x + jitter_mag/nit*2*(0.5-rand(size(x))).*x;
     neval = neval + dval;
     if  Trace
           X(:,k) = x(:);   perf(1:2,k) = [f; ng];
           perf(3:5,k-1) = [alfa; dot(h,g); dval]; end
     f_hist(k) = full(f);
     %  Check stopping criteria
     if      ng <= opts(4)
         found = 1;
     elseif  alfa*nh <= opts(5)*(opts(5) + norm(x)) 
         found = 2;
     elseif  neval >= kmax
       found = 3;
     elseif nMaxIter>0  & nit > nMaxIter
       found = 4;
%       elseif k>3 & (f_hist(1,k-1)-f_hist(1,k))/f_hist(k-1)<1e-9
%	 fprintf('Rat=%g\n',(f_hist(1,k-1)-f_hist(1,k))/f_hist(k-1));
%         found = 4;
     end
     if verbose
%    fprintf('%g>%g  %g>%g. norm(x)=%g norm(h)=%g nh=%g alfa=%g\n',ng,opts(4),alfa*nh,opts(5)*(opts(5) + norm(x)),norm(x),norm(h),nh,alfa);
    fprintf('It=%d Val=%g Alpha=%g NEval=%d Time=%g\n',nit,full(f),alfa,neval,toc);
     end
     if isnan(alfa)
       save conjgrad_drop
       keyboard
     end
   end  % iteration

%   fprintf('It=%d Val=%g NEval=%d\n',nit,full(f),neval)
   %  Set return values
   if  Trace
     X = X(:,1:k);   perf = perf(:,1:k);
   else,  X = x;  end
   info = [f  ng  alfa*nh  k-1  neval  found];

% ==========  auxiliary function  =================================

function  [x,n, f,g, opts] = check(fun,par,x0,opts0)
%  Check function call
   sx = size(x0);   n = max(sx);
   if  (min(sx) > 1)
       error('x0  should be a vector'), end
   x = x0;   [f g] = feval(fun,x,par);
   sf = size(f);   sg = size(g);
   if  any(sf-1) | ~isreal(f)
       error('f  should be a real valued scalar'), end
   if  any(sg - sx)
       error('g  should be a vector of the same type as  x'), end

   so = size(opts0);
   if  (min(so) ~= 1) | (max(so) < 9) | any(~isreal(opts0(1:9)))
       error('opts  should be a real valued vector of length 9'), end
   opts = opts0(1:9);   opts = opts(:).';
   i = find(opts(1:2) > 2);
   if  length(i)    % Set default values
     d = [2  2];   opts(i) = d(i);
   end
   i = find(opts(1:6) <= 0);
   if  length(i)    % Set default values
     d = [2  2  1  1e-4*norm(g, inf)  1e-6  5000];
%     d = [2  2  1  eps*norm(g, inf)  eps  5000];
     opts(i) = d(i);
   end
   i = find(opts(7:9) <= 0);
   if  length(i)    % Set default values
     d = [1e-6  1e-6  10; 1e-1  1e-2  5];
     %d = [1e-7  1e-7  10; 1e-5  1e-5  100];
     opts(6+i) = d(opts(2),i);
   end
return   

