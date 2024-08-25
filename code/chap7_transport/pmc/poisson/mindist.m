% Find the minimum distance from any boundary
function r=mindist(x,y,x_low,x_up,y_low,y_up)
    ymin=min((y_up-y),(y-y_low));
    xmin=min((x_up-x),(x-x_low));
    r=min(xmin,ymin);
end