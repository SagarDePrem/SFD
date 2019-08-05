function fx = ss3func1 (x)

% delta-raan at ascending node objective function

% required by sunsync3.m

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global mu ssdraan oev tnorbits raan_dot

global raan0 tiwrk yi yfinal tnode

% truncation error tolerance

tetol = 1.0e-10;

% number of differential equations

neq = 6;

% current value of orbital inclination

oev(3) = x;

% determine initial eci state vector

[rwrk, vwrk] = orb2eci(mu, oev);

% save state vector as "old" values
    
yold(1:3) = rwrk;
    
yold(4:6) = vwrk;

% propagation step size guess (seconds)

h = 60.0;

% propagation "search" step size (seconds)

dtstep = 120.0;

% initialize 'final' time

tf = 0.0;

% set initial value of objective function

% since we are starting at the ascending node
% the z component of the position vector = 0

fxnew = 0.0;

% numerically integrate and find time of nodal crossing

norbits = 0;

while(1)
    
   % initialize 'root found' indicator

   rflg = 0;

   % save current value of objective function

   fxold = fxnew;

   % set initial time to final time

   ti = tf;

   % increment final time

   tf = ti + dtstep;

   % save initial time as left side of bracket

   tisaved = ti;

   % integrate equations of motion

   ynew = rkf78('ss3eqm', neq, ti, tf, h, tetol, yold);

   % save current value of objective function
   % as the z component of the position vector

   fxnew = ynew(3);
   
   yold = ynew;
   
   % check to see if a nodal crossing
   % has been bracketed during this step

   if (fxnew * fxold < 0)
       
      rflg = 1;
      
   end   

   % if bracketed and this is an ascending node (vz > 0),
   % find time of ascending node crossing

   if (rflg == 1 && ynew(6) > 0.0)

      % load 'working' time and state vector array
      % as values on right side of bracket

      tiwrk = tf;

      yi = ynew;
      
      % root-finding convergence criterion
      
      rtol = 1.0e-8;

      % find time of ascending node crossing (troot)

      [troot, froot] = brent('ss3func2', tisaved, tiwrk, rtol);

      norbits = norbits + 1;
      
      r = yfinal(1:3);
      
      % calculate current raan (radians)
      
      ur = r / norm(r);

      raan = atan3(ur(2), ur(1));
      
      % check to see if we've integrated enough "nodal" orbits
      
      if (norbits == tnorbits)
          
         break;
         
      end
      
   end 
   
end

% "average" nodal period (seconds)

tnode = troot / norbits;

% compute objective function as difference between
% "average" regression rate and desired sun-synchronous rate

raan_dot = (raan - raan0) / troot;

fx = ssdraan - raan_dot;


