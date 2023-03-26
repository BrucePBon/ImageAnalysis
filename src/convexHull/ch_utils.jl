
function computePointsCenter( points )
    center = zeros( length(points[1]) ); 
    for p in points
        center .+= p
    end
    return Tuple( center ./ length(points) )
end

centerPoints( points, pivot ) = centerPoints!( copy(points), pivot ); 

function centerPoints!( points, pivot )

    @inbounds for idx in 1:length( points )
        vec  = points[ idx ] .- pivot; 
        dist = sqrt( sum( vec .^ 2 ) ); 
        points[idx] = vec ./ dist^2
    end
    return points
end

function convexSumTuples( vals, weights )
    convex_w = weights ./ sum(weights); 
    res = zeros( length(vals[1]) ); 
    for idx in 1:length(vals)
        @inbounds res .+= vals[idx] .* convex_w[idx]
    end
    return Tuple( res )
end
