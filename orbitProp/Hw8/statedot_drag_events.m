function [value,isterminal,direction] = statedot_drag_events(t,state)
% Locate the time when height passes through zero in a 
% decreasing direction and stop integration.

r = state(1:3);
DU = 6378.137;  %km

value = norm(r)-DU;     % Detect altitude = 0
isterminal = 1;   % Stop the integration
direction = -1;   % Negative direction only
end