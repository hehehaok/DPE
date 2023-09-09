function coorInENU = ECEF2ENU(refState, curState)
    [lat, lon, ~] = wgsxyz2lla(refState);

    slat = sind(lat);
    clat = cosd(lat);
    slon = sind(lon);
    clon = cosd(lon);
    
    R_ECEF2ENU = [ -slon,       clon,      0.0;
                   -slat*clon, -slat*slon, clat;
                    clat*clon,  clat*slon, slat];
    
    diffState = curState - refState;
    
    coorInENU = R_ECEF2ENU * diffState;
end