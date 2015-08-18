function [rate, duration] = normFiringRate(XT, YT, X_Run_Lap, Y_Run_Lap, int_at_maze,...
                               segments, show, connectgrids, roiDims,...
                               leftT, rightT, ctrNeuron, trial)
%NORMFIRINGRATE Produced a normalized unsmoothed firing rate by couting the
%               number of spikes on N segments (spatial bins) and
%               normalized by the time the animal spent in that bin.

[numLaps, N]        = size(X_Run_Lap);
rate                = zeros(N, segments, numLaps); 
duration            = zeros(segments, numLaps); 
gridsL              = get_grids(leftT, segments, connectgrids, show, roiDims);
gridsR              = get_grids(rightT, segments, connectgrids, show, roiDims);

for ilap = 1 : numLaps
   xt   = XT(int_at_maze(ilap,1):int_at_maze(ilap,2));
   yt   = YT(int_at_maze(ilap,1):int_at_maze(ilap,2));
   for icell = 1 : N
       show = 0;
       if icell == ctrNeuron; show = 1; figure(ilap), hold on;end
       
       x    = X_Run_Lap{ilap,icell};
       y    = Y_Run_Lap{ilap,icell};
      
       if strcmp(trial{ilap}, 'right') || strcmp(trial{ilap}, 'errorLeft')
           [cnt, d]= countROIspks(x, y, xt, yt, gridsR, show);
       else
           [cnt, d] = countROIspks(x, y, xt, yt, gridsL, show); 
       end
       %time normalization
       rate(icell, :, ilap) = cnt./(d);
       if icell == ctrNeuron;
            fprintf('Counting => Lap %d, cell %d, norm. f. rate %3.3f \n'...
                    ,ilap, icell, sum(rate(icell, :, ilap)))
       end
   end
   duration(:, ilap) = d;
end

end

function [count, duration]= countROIspks(x, y, xt, yt, grid, show)
    %COUNTROISPKS Counts the spikes inside a grid partition
    %   
    %   COUNT = countROIspks(x, y, grid, show)
    %           x and y are the coordinates of the spike event
    %           grid with dimension 2 x 5 x S, containes the
    %           partitions of the area in S segments. Show enables the
    %           ploting of the countign and ROI regions along with the spks
    %           xt, yt are the spatial coordinates of the animal
    %           to calculate the time spent on each
    %           bin.
    %
    %Rube Pinzon
    %version 1.0 2015
    
    xcopy   = x;
    ycopy   = y;
    numROIs = size(grid,3);
    count   = zeros(1, numROIs);
    duration    = zeros(1, numROIs);
    for iroi = 1 : numROIs
        ROI = grid(:,:,iroi);
        insideIndex = inpolygon(x,y,ROI(1,:),ROI(2, :));
        %to avoid the indefinite value, add +1 to all the bins
        duration(iroi) = sum(inpolygon(xt,yt,ROI(1,:),ROI(2, :)))+1;
        
        %remove counted spikes
        x(insideIndex) = [];
        y(insideIndex) = [];
        count(iroi)    = sum(insideIndex); 
                            
        if show
            centroid = polygonCentroid(ROI(1,:),ROI(2, :));
            text(centroid(1), centroid(2),num2str(count(iroi)), 'color', 'r')
            plot(ROI(1,:),ROI(2,:), 'r'), hold on
            plot(xt,yt)
        end   
    end
    if show
        plot(xcopy, ycopy, 'x')
    end
end

function [centroid, area] = polygonCentroid(varargin)
    %POLYGONCENTROID Compute the centroid (center of mass) of a polygon
    %
    %   CENTROID = polygonCentroid(POLY)
    %   CENTROID = polygonCentroid(PTX, PTY)
    %   Computes center of mass of a polygon defined by POLY. POLY is a N-by-2
    %   array of double containing coordinates of vertices.
    %
    %   [CENTROID AREA] = polygonCentroid(POLY)
    %   Also returns the (signed) area of the polygon. 
    %
    %   Example
    %     % Draws the centroid of a paper hen
    %     x = [0 10 20  0 -10 -20 -10 -10  0];
    %     y = [0  0 10 10  20  10  10  0 -10];
    %     poly = [x' y'];
    %     centro = polygonCentroid(poly);
    %     drawPolygon(poly);
    %     hold on; axis equal;
    %     drawPoint(centro, 'bo');
    % 
    %   References
    %   algo adapted from P. Bourke web page
    %
    %   See also:
    %   polygons2d, polygonArea, drawPolygon
    %

    %   ---------
    %   author : David Legland 
    %   INRA - TPV URPOI - BIA IMASTE
    %   created the 05/05/2004.
    %

    % Algorithme P. Bourke, vectorized version

    % HISTORY
    % 2012.02.24 vectorize code


    % parse input arguments
    if nargin == 1
        var = varargin{1};
        px = var(:,1);
        py = var(:,2);
    elseif nargin == 2
        px = varargin{1};
        py = varargin{2};
    end

    % vertex indices
    N = length(px);
    iNext = [2:N 1];

    % compute cross products
    common = px .* py(iNext) - px(iNext) .* py;
    sx = sum((px + px(iNext)) .* common);
    sy = sum((py + py(iNext)) .* common);

    % area and centroid
    area = sum(common) / 2;
    centroid = [sx sy] / 6 / area;
end

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
   
