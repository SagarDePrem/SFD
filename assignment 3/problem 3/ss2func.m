function fx = ss2func(x)

% delta-raan rate objective function

% required by sunsync2.m

% input

%  x = current orbital inclination (radians)

% output

%  fx = objective function

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global j2 j4 j22 mm drdt0 

global req2 req4 slr2 slr4 ecc2

sinc = sin(x);
   
cinc = cos(x);

% Kozai j2 + j4 perturbed mean motion (radians/second)

pmm1 = 1 + 1.5 * j2 * req2 * sqrt(1 - ecc2) * (1 - 1.5 * sinc ^ 2) / slr2;

pmm2 = (3 / 128) * j22 * req4 * sqrt(1 - ecc2) / slr4;

pmm3 = 16 * sqrt(1 - ecc2) + 25 * (1 - ecc2) - 15;

pmm4 = (30 - 96 * sqrt(1 - ecc2) - 90 * (1 - ecc2)) * cinc ^ 2;

pmm5 = (105 + 144 * sqrt(1 - ecc2) + 25 * (1 - ecc2)) * cinc ^ 4;

pmm6 = -(45 / 128) * j4 * req4 * sqrt(1 - ecc2) * ecc2 ...
       * (3 - 30 * cinc ^ 2 + 35 * cinc ^ 4) / slr4;

pmm = mm * (pmm1 + pmm2 * (pmm3 + pmm4 + pmm5) + pmm6);

% raan time rate of change

drdt1 = -1.5 * req2 * j2 * pmm * cinc / slr2;

drdt2 = 1.5 * req2 * j2 / slr2;

drdt3 = 1.5 + ecc2 / 6 - 2 * sqrt(1 - ecc2);

drdt4 = -((5 / 3) - (5 / 24) * ecc2 - 3 * sqrt(1 - ecc2)) * sinc * sinc;

drdt5 = -(35 / 8) * req4 * j4 * mm * (1 + 1.5 * ecc2) ...
        * ((12 - 21 * sinc * sinc) / 14) * cinc / slr4;

drdt = drdt1 * (1 + drdt2 * (drdt3 + drdt4)) - drdt5;

% calculate the objective function as the difference
% between the desired and current raan rate

fx = drdt - drdt0;

