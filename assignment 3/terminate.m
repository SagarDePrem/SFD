function [lookfor stop direction] = terminate(t,y)

% This function specifies the event at which ode45 terminates.
% ------------------------------------------------------------------------
lookfor = alt - 100; % � 0 when altitude � 100 km
stop = 1; % 1 means terminate at lookfor � 0; Otherwise 0
direction = -1; % -1 means zero crossing is from above
end %terminate