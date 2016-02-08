function grids = get_grids(X, segments, connect, show, roiDims)
%   GETGRIDS creates rectagular segments along a instructive path
%
%       GRIDS = get_grids(X, segments, connect, show, roiDims)
%               X is the coordinates of the instructive path along which
%               the rectagles will be palced
%               segments in the number of rectangles to fit in th elenght
%               of the accumulated distance spanned by X
%               connect (1 or 0) is an option to connect adyacent
%               rectangles or leave them disconnected.
%               show (1 or 0) to show the rectagles and the instructive
%               path. roiDims = [width lenght] of rectangles
%
%Ruben Pinzon
%version 1.0@2015
    deltaX   = diff(X(:,1));
    deltaY   = diff(X(:,2));
    accdist  = cumsum(sqrt(deltaX.^2 + deltaY.^2));
    accdist  = accdist - accdist(1);
    bin_mm   = accdist(end)/(segments+1);
       
    

    grids    = zeros(2, 5, segments);
    border_old = 1;
    for ibin = 1 : segments
        border_new = find(diff(accdist <= ibin*bin_mm));
        plot(X([border_old, border_new],1),X([border_old, border_new],2),'o-'), hold on
        
        deltaXX = X(border_new,1) - X(border_old,1);
        deltaYY = X(border_new,2) - X(border_old,2);
        offsetX  = sum(X([border_old border_new],1))/2;
        offsetY  = sum(X([border_old border_new],2))/2;
        angle    = atan2(deltaYY,deltaXX);
        
        rotation = [cos(angle) -sin(angle) 
                  sin(angle)  cos(angle)];
        translation = repmat([offsetX; offsetY],1,5);      
        ROI    = rotation*getROI(0, 0, roiDims)' + translation;
%         ROI    = rotation*ROI' + translation;

        if ibin > 1 && connect
            ROI(:,[4, 3]) = oldROI;
        end
        grids(:, :, ibin) = ROI;
        if show 
            plot(ROI(1,:),ROI(2,:), 'r')
        end
        oldROI = ROI(:,1:2);
        border_old = border_new;
        %count spikes in the grids of the central arm
    end
end

function ROI = getROI(xo, yo, roiDims)
%GETROI creates a rectagle centered at xo yo and roiDims = [width lenght]
ROI = [     
       xo + roiDims(1)/2,  yo + roiDims(2)/2;
       xo + roiDims(1)/2,  yo - roiDims(2)/2;
       xo - roiDims(1)/2,  yo - roiDims(2)/2;
       xo - roiDims(1)/2,  yo + roiDims(2)/2;
       xo + roiDims(1)/2,  yo + roiDims(2)/2
       ];
end
   
