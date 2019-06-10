function [x,flag,relres,iter,resvec] = bcgs(A,b,tol,maxit,M1,M2,x0)
%BICGSTAB BiConjugate Gradients Stabilized Method
%   X = BICGSTAB(A,B) attempts to solve the system of linear equations A*X=B
%   for X.  The coefficient matrix A must be square and the right hand
%   side (column) vector B must have length N, where A is N-by-N.  BICGSTAB
%   will start iterating from an initial guess which by default is
%   an all zero vector of length N.  Iterates are produced until the
%   method either converges, fails, or has computed the maximum number
%   of iterations.  Convergence is achieved when an iterate X has
%   relative residual NORM(B-A*X)/NORM(B) less than or equal to the
%   tolerance of the method. 

% Check for an acceptable number of input arguments
if nargin < 2
    error('at least 2 input arguments, or what are you doing here?');
elseif nargin<3
    tol = 1E-6;
    maxit = 40;
    existM1 = 0;
    existM2 = 0;
    x0 = zeros(size(b));
elseif nargin<4
    maxit = 40;
    existM1 = 0;
    existM2 = 0;
    x0 = zeros(size(b));
elseif nargin<5
    existM1 = 0;
    existM2 = 0;
    x0 = zeros(size(b));
elseif nargin<6
    existM1 = 1;
    existM2 = 0;
    x0 = zeros(size(b));
elseif nargin< 7
    existM1 = 1;
    existM2 = 1;
    x0 = zeros(size(b));
else
    existM1 = 1;
    existM2 = 1;
end

% initial parameter checks
[m,n] = size(A);
if m ~= n
  error('Matrix must be square.');
end
[mb,nb] = size(b);
if mb==m&&nb==1
    % go on
elseif mb==m&&nb==1
    b = b.';
    x0 = x0.';
else
    error('matrix and right hand side sizes do not agree!');
end
x = x0;
% Check for all zero right hand side vector => all zero solution
n2b = norm(b);                     % 2 Norm of rhs vector, b
if n2b == 0                        % if    rhs vector are all zeros
  x = zeros(n,1);                  % then  solution is all zeros
  flag = 0;                        %       a valid solution has been obtained
  relres = 0;                      %       the relative residual is 0/0
  iter = 0;                        %       no iterations need be performed
  resvec = 0;                      %       resvec(1) = norm(b-A*x) = norm(0)
  return;
end

% Set up for the method
flag = 1;
xmin = x;                          % Iterate which has minimal residual so far
imin = 0;                          % Iteration at which xmin was computed
tolb = tol * n2b;                  % Relative tolerance
r = b - A * x;                     % Zero-th residual
normr = norm(r);                   % Norm of residual

if normr <= tolb                   % Initial guess is a good enough solution
  flag = 0;
  relres = normr / n2b;
  iter = 0;
  resvec = normr;
  return;
end

rt = r;                            % Shadow residual
resvec = zeros(2*maxit+1,1);       % Preallocate vector for norm of residuals
resvec(1) = normr;                 % resvec(1) = norm(b-A*x0)
normrmin = normr;                  % Norm of residual from xmin
rho = 1;
omega = 1;
stag = 0;                          % stagnation of the method

% loop over maxit iterations (unless convergence or failure)

for i = 1 : maxit
  rho1 = rho;
  rho = rt' * r;
  if rho == 0.0 || rho == Inf
    flag = 4;
    resvec = resvec(1:2*i-1);
    break;
  end
  if i == 1
    p = r;
  else
    beta = (rho/rho1)*(alph/omega);
    if beta == 0 || beta == Inf || isnan(beta)
      flag = 4;
      break;
    end
    p = r + beta * (p - omega * v);
  end
  if existM1==1
    ph1 = M1 \ p;
    if isinf(norm(ph1,inf))
      flag = 2;
      resvec = resvec(1:2*i-1);
      break;
    end
  else
    ph1 = p;
  end
  if existM2==1
    ph = M2 \ ph1;
    if isinf(norm(ph,inf))
      flag = 2;
      resvec = resvec(1:2*i-1);
      break;
    end
  else
    ph = ph1;
  end
  v = A * ph;
  rtv = rt' * v;
  if rtv == 0 || rtv == Inf
    flag = 4;
    resvec = resvec(1:2*i-1);
    break;
  end
  alph = rho / rtv;
  if alph == Inf
    flag = 4;
    resvec = resvec(1:2*i-1);
    break;
  end
  if alph == 0                    % stagnation of the method
    stag = 1;
  end

% Check for stagnation of the method
  if stag == 0
    stagtest = zeros(n,1);
    ind = (x ~= 0);
    stagtest(ind) = ph(ind) ./ x(ind);
    stagtest(~ind & ph ~= 0) = Inf;
    if abs(alph)*norm(stagtest,inf) < eps
      stag = 1;
    end
  end

  xhalf = x + alph * ph;          % form the "half" iterate
  normr = norm(b-A*xhalf);
  resvec(2*i) = normr;

  if normr <= tolb                 % check for convergence
    x = xhalf;
    flag = 0;
    iter = i - 0.5;
    resvec = resvec(1:2*i);
    break;
  end

  if stag == 1
    flag = 3;
    resvec = resvec(1:2*i);
    break;
  end

  if normr < normrmin              % update minimal norm quantities
    normrmin = normr;
    xmin = xhalf;
    imin = i - 0.5;
  end

  s = r - alph * v;               % residual associated with xhalf
  if existM1==1
    sh1 = M1 \ s;
    if isinf(norm(sh1,inf))
      flag = 2;
      resvec = resvec(1:2*i);
      break;
    end
  else
    sh1 = s;
  end
  if existM2==1
    sh = M2 \ sh1;
    if isinf(norm(sh,inf))
      flag = 2;
      resvec = resvec(1:2*i);
      break;
    end
  else
    sh = sh1;
  end
  t = A * sh;
  tt = t' * t;
  if tt == 0 || tt == Inf
    flag = 4;
    resvec = resvec(1:2*i);
    break;
  end
  omega = (t' * s) / tt;
  if omega == Inf
    flag = 4;
    resvec = resvec(1:2*i);
    break;
  end
  if omega == 0                    % stagnation of the method
    stag = 1;
  end

% Check for stagnation of the method
  if stag == 0
    stagtest = zeros(n,1);
    ind = (xhalf ~= 0);
    stagtest(ind) = sh(ind) ./ xhalf(ind);
    stagtest(~ind & sh ~= 0) = Inf;
    if abs(omega)*norm(stagtest,inf) < eps
      stag = 1;
    end
  end

  x = xhalf + omega * sh;          % x = (x + alpha * ph) + omega * sh
  normr = norm(b-A*x);
  resvec(2*i+1) = normr;

  if normr <= tolb                 % check for convergence
    flag = 0;
    iter = i;
    resvec = resvec(1:2*i+1);
    break;
  end

  if stag == 1
    flag = 3;
    resvec = resvec(1:2*i+1);
    break;
  end

  if normr < normrmin              % update minimal norm quantities
    normrmin = normr;
    xmin = x;
    imin = i;
  end

  r = s - omega * t;

end                                % for i = 1 : maxit

% returned solution is first with minimal residual
%if flag == 0
  relres = normr / n2b;
  iter = i;
%else
%   x = xmin;
%   iter = imin;
%   relres = normrmin / n2b;
% end

return