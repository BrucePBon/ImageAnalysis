using PyPlot

function convexHull_reverse( points, center )
    centered_points = centerPoints( points, center )
    edges, norms    = convexHull( centered_points ); 
    return edges, norms
end

function convexHull_reverse_mask( points, center; repval=false )
    centered_points = centerPoints( points, center )
    edges, norms, mask = convexHull_mask( centered_points, repval=repval ); 
    return edges, norms, mask
end

function convexHull_reverse_iterative( points, pivot; c=nothing )

    # Transforming points to compute the reversed convex hull around the pivot
    trans_points = centerPoints( points, pivot ); 
    
    # Preparing arrays for computing the initial reversed hull
    Npoints = length(trans_points); 
    mask    = zeros( Bool, Npoints ); 
    edges   = Array{Tuple{Int64,Int64},1}(undef,Npoints*3);
    norms   = Array{Tuple{Float64,Float64},1}(undef,Npoints*3);
    visible = Array{Int64,1}(undef,Npoints*2); # visibleEdges

    # Computing the initial reversed hull
    lastEdge, center = initiateHullMask!( trans_points, edges, norms, mask );
    lastEdge = expandHullMask!( lastEdge, center, trans_points, edges, norms, visible, mask, repval=false )

    # For each edge on the hull, compute a new pivot on that edge, and expand the reverse hull
    for iters in 1:4

        edge = edges[ rem( iters, lastEdge) + 1 ];

        # New pivot and recentering points
        p1, p2       = points[edge[1]], points[edge[2]]; 
        new_pivot    = convexSumTuples( ( p1, p2, pivot ), (1, 1, 2) ); 
        trans_points = centerPoints( points, new_pivot )
        tp1, tp2     = trans_points[edge[1]], trans_points[edge[2]]; 

        
        # debuging
        debug = true
        if debug
            figure( nothing, (10, 4) )
            title( string(iters) )
            subplot( 1, 4, 1 ); 
            ax = PyPlot.gca();     
            ax[:set_aspect]("equal")
            scatter( [ p[2] for p in points ], [ p[1] for p in points ], c=c ); 
            for edge in edges[1:lastEdge]
                plot( [ points[edge[1]][2], points[edge[2]][2] ], [ points[edge[1]][1], points[edge[2]][1] ], "k-" )
            end
            scatter( [ p1[2], p2[2] ], [ p1[1], p2[1] ] ); 
            scatter( [ new_pivot[2] ], [ new_pivot[1] ] ); 

            subplot( 1, 4, 2 ); 
            ax = PyPlot.gca();     
            ax[:set_aspect]("equal")
            scatter( [ p[2] for p in trans_points ], [ p[1] for p in trans_points ], c=c )
            scatter( [ tp1[2], tp2[2] ], [ tp1[1], tp2[1] ] ); 

            for edge in edges[1:lastEdge]
                plot( [ trans_points[edge[1]][2], trans_points[edge[2]][2] ], [ trans_points[edge[1]][1], trans_points[edge[2]][1] ], "k-" )
            end
        end

        # recomputing normals
        for eidx in 1:lastEdge
            norms[eidx] = edgeNormal( trans_points, edges[eidx], (0,0) );
        end

        lastEdge = expandHullMask!( lastEdge, center, trans_points, edges, norms, visible, mask, repval=true ) 

        debug = true
        if debug
            subplot( 1, 4, 3 ); 
            ax = PyPlot.gca();     
            ax[:set_aspect]("equal")
            scatter( [ p[2] for p in trans_points ], [ p[1] for p in trans_points ], c=c )
            scatter( [ tp1[2], tp2[2] ], [ tp1[1], tp2[1] ] ); 

            for edge in edges[1:lastEdge]
                plot( [ trans_points[edge[1]][2], trans_points[edge[2]][2] ], [ trans_points[edge[1]][1], trans_points[edge[2]][1] ], "k-" )
            end

            subplot( 1, 4, 4 ); 
            ax = PyPlot.gca();     
            ax[:set_aspect]("equal")
            scatter( [ p[2] for p in points ], [ p[1] for p in points ], c=c ); 
            for edge in edges[1:lastEdge]
                plot( [ points[edge[1]][2], points[edge[2]][2] ], [ points[edge[1]][1], points[edge[2]][1] ], "k-" )
            end
            scatter( [ p1[2], p2[2] ], [ p1[1], p2[1] ] ); 
            scatter( [ new_pivot[2] ], [ new_pivot[1] ] ); 


        end
    end

    return edges[1:lastEdge], norms[1:lastEdge]
end


#= ITERATIVE REVERSE HULL PSEUDOCODE

    edges, normals, mask = compute reverse hull with mask

    points in hull = sum( mask ); 

    points added to hull = points in hull; 

    # ITERATE UNTIL NO MORE POINTS ARE ADDED TO THE HULL
    while points added to hull > 0
        
        # FOR EACH EDGE OF THE HULL, SHIFT THE CENTER TO THAT EDGE, RECOMPUTE COORDINATES AND normals
        # AND TRY TO EXPAND THE HULL
        for each edge in the hull
        
            center = edge midpoint

            "reverse"" point coordinates w/resp to center and recompute convex hull normals

            edges, normals, mask = expand reverse hull with mask

        end

        points added to hull = sum( mask ) - points in hull; 

        points in hull = sum( mask ); 
    end

=#