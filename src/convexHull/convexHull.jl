""" 2D Implementation """
# Look at the 3D implementation for commented code. 

function convexHull( points::Array{Tuple{T,T},1} ) where {T<:Real}

    @assert length(points) >= 3

    edges = Array{Tuple{Int64,Int64},1}(undef,length(points)*3);
    norms = Array{Tuple{Float64,Float64},1}(undef,length(points)*3);
	visibleEdges = Array{Int64,1}(undef,length(points)*2);

	lastEdge, center, _, _ = initiateHull!( points, edges, norms )

	lastEdge, _, _ = expandHull!( lastEdge, center, points, edges, norms, visibleEdges )
	
    return edges[1:lastEdge], norms[1:lastEdge]
end

##### Functions for 2D ConvexHull

function initiateHull!( points::Array{Tuple{T,T},1}, edges, norms ) where {T}
    
	prepareStartingPoints!( points )

    edges[1:3] = [ ( 1, 2 ), (2, 3), ( 3, 1 ) ];
    for eidx in 1:3 
		norms[eidx] = edgeNormal( points, edges[eidx], sum( points[1:3] ) ./ 3 )
    end

	return 3, sum( points[1:3] ) ./ 3, edges, norms
end

function prepareStartingPoints!( points::Array{Tuple{T,T},1} ) where {T}

	idx = 2
	while areEqual( points[1], points[idx] ) && idx < length(points)
		idx += 1;
	end
	swap_points!( points, 2, idx ); 

	idx = 3
	while areColinear( points[1], points[2], points[idx] ) && idx < length(points)
		idx += 1;
	end
	swap_points!( points, 3, idx ); 
end

function edgeNormal( points, edge, center )
	p1, p2 = points[edge[1]], points[edge[2]]; 
	norm   = normcross_( p1, p2 ); 
	vec    = center .- p1;
	return sign( dot_( vec, norm ) ) .* norm;
end

function findVisibleEdges!( lastEdge, test_point, points, edges, norms, visibleEdges )

	numVisibleEdges = 0;
	# Going through each face
	for eidx in 1:lastEdge;
		# vector from a point in the current face to new_point
		edge2point_vec = test_point .- points[edges[eidx][1]];
		edge_normal    = norms[ eidx ]; 
		v_normal_dot   = dot_( edge2point_vec, edge_normal )
		# if "v_normal_dot" is positive, add the face index to "visibleFaces". 
		if sign( v_normal_dot ) == 1
			numVisibleEdges += 1;
			visibleEdges[numVisibleEdges] = eidx;
		end
	end
	return numVisibleEdges
end

function findRepeatedPoints( numVisibleEdges, visibleEdges, edges )
	boundary    = Array{Int64,1}(undef, numVisibleEdges*2);
	isRepeated  = zeros( Bool, numVisibleEdges*3 );
	pointidx = 1
	@inbounds for edge in edges[ visibleEdges[1:numVisibleEdges] ]
		# for each visible edge --> go through its 2 points
		for p in [ edge[1], edge[2] ]
			# for each of its 2 points --> store it to "boundary" and check if it is not repeated
			boundary[pointidx] = p
			for prev_pointidx in 1:pointidx-1
				# each previous boundary --> check that the current point has not already been added
				prev_point = boundary[prev_pointidx]
				if p == prev_point
					isRepeated[  pointidx   ] = true;
					isRepeated[prev_pointidx] = true;
				end
			end
			pointidx += 1 
		end 
	end 
	return boundary, isRepeated
end

function addNewEdges!( points, edges, norms, visibleEdges, newpoint_idx, lastEdge, center, numVisibleEdges )

	# Isolating non-repeated edges among visibleFace.
	boundary, isRepeated = findRepeatedPoints( numVisibleEdges, visibleEdges, edges ); 
	# For each boundary point
	counter = 0;
	for bidx in 1:length(boundary)
		# For each non-repeated edge -> create a new face between its points and "new_point"
		if isRepeated[bidx]
			continue
		else
			counter += 1
			# Instead of simply adding the newfaces at the end, we place the first $(numvisibleFaces) newfaces
			# in the position of the visibleFaces, as these have to be deleted from the convex hull anyways. If 
			# numVisibleFaces is exhausted (counter > numVisibleFaces), we add at the end. 
			newpointidx = ( counter <= numVisibleEdges ) ? visibleEdges[counter] : lastEdge + (counter-numVisibleEdges);
			# Adding new face and its normal vector
			edges[newpointidx] = ( boundary[bidx], newpoint_idx );
			norms[newpointidx] = edgeNormal( points, edges[newpointidx], center )
		end
	end

	remaining = numVisibleEdges - counter; 

	# "remaining" == 0
	#     All visible faces have been replaced by all new faces. "lastFace" stays the same
	# "remaining" < 0
	#     There were more newfaces than visibleFaces, and "lastFace" has to increase. 
	# "remaining" > 0
	#     There are visibleFaces that have not been replaced and need to be removed. "lastFaces" has to decrese.
	#     The loop below only runs if (remaining > 0), so (numVisibleFaces-remaining+1) <= numVisibleFaces

	offset = 0;
	for visible_eidx in visibleEdges[ numVisibleEdges-remaining+1:numVisibleEdges ]
		# For every remaining face that we remove, all subsquent indices are offset by 1. "offset" accounts for that.
		for eidx in visible_eidx-offset:lastEdge-offset
			edges[eidx] = edges[eidx+1]
			norms[eidx] = norms[eidx+1]
		end
		offset += 1;
	end

	return lastEdge - remaining
end

function expandHull!( lastEdge, center, points::Array{Tuple{T,T},1}, edges, norms, visibleEdges ) where {T}
	
	@inbounds for pidx in 4:length(points)
		new_point = points[pidx];
		numVisibleEdges = findVisibleEdges!( lastEdge, new_point, points, edges, norms, visibleEdges )
		lastEdge = addNewEdges!( points, edges, norms, visibleEdges, pidx, lastEdge, center, numVisibleEdges );
	end 
	return lastEdge, edges, norms
end