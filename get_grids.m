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

    deltaX   = diff(X(1:segments:end,1));
    deltaY   = diff(X(1:segments:end,2));
    X        = X(1:segments:end,:);
    angle    = atan2(deltaY,deltaX);

    grids    = zeros(2, 5, length(deltaX)-1);
    for ibin = 1 : length(deltaX)-1

        plot(X(ibin:ibin+1,1),X(ibin:ibin+1,2),'o-'), hold on

        rotation = [cos(angle(ibin)) -sin(angle(ibin)) 
                  sin(angle(ibin))  cos(angle(ibin))];
        translation = repmat([sum(X(ibin:ibin+1,1))/2; sum(X(ibin:ibin+1,2))/2],1,5);      
        ROI    = rotation*getROI(0, 0, roiDims)' + translation;
        if ibin > 1 && connect
            ROI(:,[4, 3]) = oldROI;
        end
        grids(:, :, ibin) = ROI;
        if show 
            plot(ROI(1,:),ROI(2,:), 'r')
        end
        oldROI = ROI(:,1:2);
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
   
