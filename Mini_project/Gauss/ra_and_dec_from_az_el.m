% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
function [ra,dec] = ra_and_dec_from_az_el(az,el,lat,theta) 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%{   
  This function calculates the right ascension and the 
  declination from the geocentric equatorial position vector 
  
  r       - position vector 
  l, m, n - direction cosines of r 
  ra      - right ascension (degrees) 
  dec     - declination (degrees) 
%} 
% ---------------------------------------------- 
dec = asind(cosd(lat)*cosd(az)*cosd(el)+sind(lat)*sind(el)); 
  
if az > 0 && az < 180 
    h = 360 - acosd((cosd(lat)*sind(el)-sind(lat)*cosd(az)*cosd(el))/cosd(dec)); 
else 
    h = acosd((cosd(lat)*sind(el)-sind(lat)*cosd(az)*cosd(el))/cosd(dec)); 
end 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
ra = theta - h; 