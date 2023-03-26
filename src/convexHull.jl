# Transforming DxN point arrays into arrays of tuples, for performance reasons. 
function convexHull( points::Array{T,2} ) where {T<:Real}
	if     size(points,1) == 2
		return convexHull2D( DN2Tuple(points) ); 
	elseif size(points,2) == 3
		return convexHull3D( DN2Tuple(points) ); 
	end
end

""" 2D Implementation """

function convexHull( points::Array{Tuple{T,T},1} ) where {T<:Real}

	Npoints = length(points); 
    @assert Npoints >= 3

	# preallocating memory to store necessary arrays for storing visibleEdges, hullEdges and hullNormals.
	visible = Array{Int64,1}(undef, Npoints*2);
    edges   = Array{Tuple{Int64,Int64},1}(undef, Npoints*3);
    norms   = Array{Tuple{Float64,Float64},1}(undef, Npoints*3);

	# Sampling the initial triangle and expanding the convex hull. 
	lastEdge, center = initiateHull!( points, edges, norms )
	lastEdge         = expandHull!( lastEdge, center, points, edges, norms, visible )
	
    return edges[1:lastEdge], norms[1:lastEdge]
end

##### Functions for 2D ConvexHull

function initiateHull!( points::Array{Tuple{T,T},1}, edges, norms ) where {T}

	center = prepareNonCollinearPoints!( points ); 

	edges[1:3] .= [ ( 1, 2 ), ( 2, 3 ), ( 3, 1 ) ];
    for eidx in 1:3
		norms[eidx] = edgeNormal( points, edges[eidx], center );
    end
	return 3, center
end

function prepareNonCollinearPoints!( points::Array{Tuple{T,T},1} ) where {T}

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

	( areEqual(points[1:2]...) || areColinear(points[1:3]...) ) && error("All input points are equal or colinear.")
	
	return computePointsCenter( points[1:3] ); 
end

function edgeNormal( points, edge, center )

	p1, p2 = points[edge[1]], points[edge[2]]; 
	norm   = normcross_( p1, p2 ); 
	vec    = center .- p1;
	return -1 * sign(dot_(vec,norm)) .* norm; #( sign( dot_( vec, norm ) ) == -1 ) ? norm : -1 .* norm;
end

function findVisibleEdges!( lastEdge, pidx, points, edges, norms, visibleEdges )

	numVisibleEdges = 0;
	# Going through each face
	for eidx in 1:lastEdge;
		# vector from P1 the current edge to new_point
		edge2point_vec = points[pidx] .- points[edges[eidx][1]];
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
		numVisibleEdges = findVisibleEdges!( lastEdge, pidx, points, edges, norms, visibleEdges )
		lastEdge = addNewEdges!( points, edges, norms, visibleEdges, pidx, lastEdge, center, numVisibleEdges );
	end 
	
	return lastEdge
end

# Convex hull algorithm accepting a mask of points that "have already been processed" or "need to be processed". this
# adaptation exists in order to make iterative reverse convex Hulls more practical. 

function convexHullMask( points::Array{Tuple{T,T},1}; repval=true ) where {T<:Real}

	Npoints = length( points ); 
    @assert ( Npoints >= 3 )

	mask    = zeros( Bool, Npoints ); 
    edges   = Array{Tuple{Int64,Int64},1}(undef, Npoints*3);
    norms   = Array{Tuple{Float64,Float64},1}(undef, Npoints*3);
	visible = Array{Int64,1}(undef, Npoints*2);

	lastEdge, center = initiateHullMask!( points, edges, norms, mask )
	lastEdge = expandHullMask!( lastEdge, center, points, edges, norms, visible, mask, repval=repval )
	
    return edges[1:lastEdge], norms[1:lastEdge], mask
end

function initiateHullMask!( points, edges, norms, mask )
	mask[1:3] .= true
	return initiateHull!( points, edges, norms )
end

function expandHullMask!( lastEdge, center, points::Array{Tuple{T,T},1}, edges, norms, visibleEdges, mask; repval=true ) where {T}
	
	for pidx in 4:length(points)
		( mask[pidx] ) && ( continue; )
		numvisible = findVisibleEdges!( lastEdge, pidx, points, edges, norms, visibleEdges )
		lastEdge   = addNewEdgesMask!( points, edges, norms, visibleEdges, mask, pidx, lastEdge, center, numvisible, repval=repval );
	end 
	return lastEdge
end

function addNewEdgesMask!( points, edges, norms, visibleEdges, mask, point_idx, lastEdge, center, numVisibleEdges; repval=false )

	# Isolating non-repeated edges among visibleFace.
	boundary, isRepeated = findRepeatedPoints( numVisibleEdges, visibleEdges, edges ); 
	# For each boundary point
	counter = 0;
	for bidx in 1:length(boundary)
		# For each non-repeated edge -> create a new face between its points and "new_point"
		if isRepeated[bidx]
			mask[boundary[bidx]] = repval
			# If a point has been taken out of the hull, you don't want to process it again...
			continue
		else
			counter += 1
			# Instead of simply adding the newfaces at the end, we place the first $(numvisibleFaces) newfaces
			# in the position of the visibleFaces, as these have to be deleted from the convex hull anyways. If 
			# numVisibleFaces is exhausted (counter > numVisibleFaces), we add at the end. 
			newpointidx = ( counter <= numVisibleEdges ) ? visibleEdges[counter] : lastEdge + (counter-numVisibleEdges);
			# Adding new face and its normal vector
			edges[newpointidx] = ( boundary[bidx], point_idx );
			norms[newpointidx] = edgeNormal( points, edges[newpointidx], center )
			   mask[point_idx] = true
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


""" Animated version of 2D Convex Hull """

function convexHull_animate( points::Array{Tuple{T,T},1} ) where {T<:Real}

    @assert length(points) >= 3
	
	edge_animation = []

    edges = Array{Tuple{Int64,Int64},1}(undef,length(points)*3);
    norms = Array{Tuple{Float64,Float64},1}(undef,length(points)*3);
	visibleEdges = Array{Int64,1}(undef,length(points)*2);

	prepareStartingPoints!( points )
    edges[1] = ( 1, 2 );
    edges[2] = ( 2, 3 );
    edges[3] = ( 3, 1 );
    lastEdge = 3;

	center = ( points[1] .+ points[2] .+ points[3] ) ./ 3;
    for eidx in 1:3
		norms[eidx] = edgeNormal( points, edges[eidx], center )
    end

	expandHull_anim!( lastEdge, center, points, edges, norms, visibleEdges, edge_animation ); 
	
    return edge_animation
end

function expandHull_anim!( lastEdge, center, points::Array{Tuple{T,T},1}, edges, norms, visibleEdges, edge_anim ) where {T}
	
	@inbounds for pidx in 4:length(points)
		new_point = points[pidx];
		numVisibleEdges = findVisibleEdges!( lastEdge, new_point, points, edges, norms, visibleEdges )
		lastEdge = addNewEdges!( points, edges, norms, visibleEdges, pidx, lastEdge, center, numVisibleEdges );
		push!( edge_anim, edges[1:lastEdge] ); 
	end 
	return lastEdge, edge_anim
end

### Masked convex hull animation


function convexHullMask_animate( points::Array{Tuple{T,T},1}; repval=true ) where {T<:Real}

    @assert length(points) >= 3
	
	edge_animation = []
	mask_animation = []

    edges = Array{Tuple{Int64,Int64},1}(undef,length(points)*3);
    norms = Array{Tuple{Float64,Float64},1}(undef,length(points)*3);
	visibleEdges = Array{Int64,1}(undef,length(points)*2);
	mask = zeros( Bool, length(points) );

	lastEdge, center = initiateHull!( points, edges, norms )
	mask[1:3] .= true

	expandHullMask_anim!( lastEdge, center, points, edges, norms, visibleEdges, mask, edge_animation, mask_animation, repval=repval ); 
	
    return edge_animation, mask_animation
end

function expandHullMask_anim!( lastEdge, center, points::Array{Tuple{T,T},1}, edges, norms, visibleEdges, mask, edge_anim, mask_anim; repval=true ) where {T}
	
	@inbounds for pidx in 4:length(points)
		if mask[pidx]
			continue; #already in the convex hull / don't need to be processed
		end
		new_point = points[pidx];
		numVisibleEdges = findVisibleEdges!( lastEdge, new_point, points, edges, norms, visibleEdges )
		lastEdge = addNewEdgesMask!( points, edges, norms, visibleEdges, mask, pidx, lastEdge, center, numVisibleEdges, repval=repval );
		push!( edge_anim, edges[1:lastEdge] ); 
		push!( mask_anim, copy(mask) ); 

	end 
	return lastEdge, edge_anim, mask_anim
end





""" 3D Implementation """



# Convex hull inputs accepting an array of Tuples.

function convexHull( points::Array{Tuple{T,T,T},1} ) where {T<:Real}

	Npoints = length(points); 
    @assert Npoints >= 4

    # Initializing the arrays to store "faces", "normals" and "visibleFaces" with more
	# than enough space $(length(points)*3). I do this because I believe that having a 
	# dynamic list of faces, and using push! and pop! is slower than working on fixed-size
	# arrays. TODO: test this preconception
	visible = Array{Int64,1}(undef, Npoints*4);
    faces   = Array{Tuple{Int64,Int64,Int64},1}(undef, Npoints*5);
    norms   = Array{Tuple{Float64,Float64,Float64},1}(undef, Npoints*5);

	lastFace, center = initiateHull!( points, faces, norms )

	# Going over the rest of points, and expanding the convexHull, if necessary.
	lastFace = expandHull!( lastFace, center, points, faces, norms, visible );

    return faces[1:lastFace], norms[1:lastFace]
end

function convexHull_animate( points::Array{Tuple{T,T,T},1} ) where {T<:Real}

	Npoints = length(points); 
    @assert Npoints >= 4

	face_animation = []
	norm_animation = []

	visible = Array{Int64,1}(undef, Npoints*5);
    faces   = Array{Tuple{Int64,Int64,Int64},1}(undef, Npoints*5);
    norms   = Array{Tuple{Float64,Float64,Float64},1}(undef, Npoints*5);

	lastFace, center = initiateHull!( points, faces, norms )

	push!( face_animation, faces[1:lastFace] )
	push!( norm_animation, norms[1:lastFace] )

	# Going over the rest of points, and expanding the convexHull, if necessary.
	expandHull_anim!( lastFace, center, points, faces, norms, visible, face_animation, norm_animation );

    return face_animation, norm_animation
end

function initiateHull!( points::Array{Tuple{T,T,T},1}, faces, norms ) where {T}

	center = prepareNonCoplanarPoints!( points ); 

    faces[1:4] .= [ ( 1, 2, 3 ), ( 2, 3, 4 ), ( 3, 4, 1 ), ( 4, 1, 2 ) ];
    for fidx in 1:4
		norms[fidx] = faceNormal( points, faces[fidx], center );
    end
	return 4, center
end

function prepareNonCoplanarPoints!( points::Array{Tuple{T,T,T},1} ) where {T}

	idx = 2
	while areEqual( points[1], points[idx] ) && idx <= length(points)
		idx += 1;
	end
	( idx > length(points) ) && error( "All same points." )
	swap_points!( points, 2, idx ); 

	idx = 3
	while areColinear( points[1], points[2], points[idx] ) && idx <= length(points)
		idx += 1;
	end
	( idx > length(points) ) && error( "All colinear." )
	swap_points!( points, 3, idx ); 

	idx = 4

	while areCoplanar( points[1], points[2], points[3], points[idx] ) && idx <= length(points)
		idx += 1;
	end
	( idx > length(points) ) && error( "All coplanar." )
	swap_points!( points, 4, idx ); 

	return computePointsCenter( points[1:4] ); 
end

"""
	Compute the normal vector for an edge/face, ensuring that it points away from the convex hull's center.
"""
function faceNormal( points, face, center )
	p1, p2, p3 = points[face[1]], points[face[2]], points[face[3]]; 
	norm       = normcross_( p2 .- p1, p3 .- p1 ); 
	vec        = center .- p1;
	return -1 * sign(dot_(vec,norm)) .* norm;
end

""" 
	Checking the number of edges/faces of the convex hull "visible" to "new_point". 
""" 
function findVisibleFaces!( lastFace, pidx, points, faces, norms, visibleFaces )

	numVisibleFaces = 0;
	# Going through each face
	for fidx in 1:lastFace;
		# vector from a point in the current face to new_point
		face2point_vec = points[pidx] .- points[faces[fidx][1]];
		face_normal    = norms[ fidx ]; 
		v_normal_dot   = dot_( face2point_vec, face_normal )
		# if "v_normal_dot" is positive, add the face index to "visibleFaces". 
		if sign( v_normal_dot ) == 1
			numVisibleFaces += 1;
			visibleFaces[numVisibleFaces] = fidx;
		end
	end
	return numVisibleFaces
end

"""
	Among the visible faces, collecting their edges and marking those edges which are
	repeated. Those edges that are not repeated form the "boundary" between visible 
	faces and non-visible faces. New faces will be built between "new_point" and the 
	points on the "boundary".
"""

function findRepeatedEdges( numVisibleFaces, visibleFaces, faces )

	boundary    = Array{Tuple{Int64,Int64},1}(undef, numVisibleFaces*3);
	isRepeated  = zeros( Bool, numVisibleFaces*3 );
	edgeidx     = 0
	@inbounds for face in faces[ visibleFaces[1:numVisibleFaces] ]
		# for each visible face --> go through its 3 edges
		for (p1, p2) in [ (face[1],face[2]), (face[2],face[3]), (face[3],face[1]) ]
			# for each of its 3 edges --> store it to "boundary" and check if it is not repeated
			edgeidx += 1 
			boundary[edgeidx] = (p1,p2)
			for prev_edgeidx in 1:edgeidx-1
				# each previous boundary --> check that the current edge has not already been added
				prev_edge = boundary[prev_edgeidx]
				if (p1,p2) == prev_edge || (p2,p1) == prev_edge
					isRepeated[  edgeidx   ] = true;
					isRepeated[prev_edgeidx] = true;
				end
			end
		end 
	end 

	return boundary, isRepeated
end


"""
	Computes visible faces and builds new faces from the newpoint and the points on the boundary of
	visible faces.
"""

function addNewFaces!( points, faces, norms, visibleFaces,
					   newpoint_idx, lastFace, center, numVisibleFaces                       
                     )

	# Isolating non-repeated edges among visibleFace.
	boundary, isRepeated = findRepeatedEdges( numVisibleFaces, visibleFaces, faces ); 

	# For each boundary edge
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
			newfaceidx = ( counter <= numVisibleFaces ) ? visibleFaces[counter] : lastFace + (counter-numVisibleFaces);
			# Adding new face and its normal vector
			faces[newfaceidx] = ( boundary[bidx][1], boundary[bidx][2], newpoint_idx );
			norms[newfaceidx] = faceNormal( points, faces[newfaceidx], center )
		end
	end

	remaining = numVisibleFaces - counter; 

	# "remaining" == 0
	#     All visible faces have been replaced by all new faces. "lastFace" stays the same
	# "remaining" < 0
	#     There were more newfaces than visibleFaces, and "lastFace" has to increase. 
	# "remaining" > 0
	#     There are visibleFaces that have not been replaced and need to be removed. "lastFaces" has to decrese.
	#     The loop below only runs if (remaining > 0), so (numVisibleFaces-remaining+1) <= numVisibleFaces
	if remaining > 0 
		offset = 0;
		for visible_fidx in visibleFaces[ numVisibleFaces-remaining+1:numVisibleFaces ]
			# For every remaining face that we remove, all subsquent indices are offset by 1. "offset" accounts for that.
			for fidx in visible_fidx-offset:lastFace-1-offset
				faces[fidx] = faces[fidx+1]
				norms[fidx] = norms[fidx+1]
			end
			offset += 1;
		end
	end

	return lastFace - remaining
end

function expandHull!( lastFace, center, points::Array{Tuple{T,T,T},1}, faces, norms, visibleFaces ) where {T}
	
	@inbounds for pidx in 5:length(points)
		numVisibleFaces = findVisibleFaces!( lastFace, pidx, points, faces, norms, visibleFaces )
		lastFace = addNewFaces!( points, faces, norms, visibleFaces, pidx, lastFace, center, numVisibleFaces );
	end 
	
	return lastFace
end


function expandHull_anim!( lastFace, center, points::Array{Tuple{T,T,T},1}, faces, norms, visibleFaces, face_anim, norm_anim ) where {T}
	
	@inbounds for pidx in 5:length(points)
		numVisibleFaces = findVisibleFaces!( lastFace, pidx, points, faces, norms, visibleFaces )
		lastFace = addNewFaces!( points, faces, norms, visibleFaces, pidx, lastFace, center, numVisibleFaces );
		push!( face_anim, faces[1:lastFace] ); 
		push!( norm_anim, norms[1:lastFace] ); 
	end 
	return lastFace, face_anim, norm_anim
end

# To start the expansion of the convex hull we need, at
# least, (2D) three non-equal and non-colinear points or 
# (3D) four non-equal, non-colinear and non-coplanar points. 

areEqual( x, y ) = ( x == y )
areColinear( x, y, z ) = abs( normdot_( x .- z, y .- z ) ) > 0.999;
function areCoplanar( x, y, z, p )
	n1 = cross_( x .- y, x .- z ); 
	n2 = cross_( x .- y, x .- p );
	return abs( normdot_( n1, n2 ) ) > 0.95
end

function swap_points!( points, idx0, idx1 )
	if ( idx0 !== idx1 )
		temp = points[idx0]
		points[idx0] = points[idx1]
		points[idx1] = temp
	end
end































function convexHullSIMP2( points::Array{Tuple{T,T,T},1}; typ=UInt16 ) where {T<:Real}

    @assert length(points) >= 4
    # length(points) for now, to have a static size

    faces = Array{Tuple{typ,typ,typ},1}(undef,0);
    norms = Array{Tuple{Float64,Float64,Float64},1}(undef,0);


    push!( faces, (1,2,3) )
	push!( faces, (2,3,4) )
	push!( faces, (3,4,1) )
	push!( faces, (4,1,2) )

	center   = ( points[1] .+ points[2] .+ points[3] .+ points[4] ) ./ 4;
    for fidx in 1:length(faces)
		push!( norms, faceNormalCenterSIMP( points, faces[fidx], center ) )
    end

	@inbounds for pidx in 5:length(points);

		new_point = points[pidx];

		visiblefaces = Array{typ,1}(undef,0)

		for fidx in 1:length(faces);

			vp = new_point .- points[faces[fidx][1]];
			if sign( sum( vp .* norms[ fidx ] ) ) == 1
				push!( visibleFaces, fidx );
			end
		end

		numVisibleFaces = length(visibleFaces);

		lastFace = addNewFacesSIMP2!( points, faces, norms, center, lastFace,
						             pidx, visibleFaces, numVisibleFaces );

	end # each input point

    return faces[1:lastFace], norms[1:lastFace]
end

"""
	needs to be tested
"""
# can this be made more compact?
@inline function addNewFacesSIMP2!( points,
                               faces::Array{Tuple{T,T,T},1},
							   norms::Array{Tuple{Float64,Float64,Float64},1},
							   center,
							   lastFace,
							   newpointindex,
							   visibleFaces::Array{T,1},
                               numVisibleFaces
                             ) where {T<:Integer}

	boundary   = Array{Tuple{T,T},1}(undef, numVisibleFaces*3);
	isRepeated = zeros( Bool, numVisibleFaces*3 );

	lidx, nrep = 1, 0;

	@inbounds for face in faces[ visibleFaces[1:numVisibleFaces] ]
		edges = [(face[1],face[2]), (face[2],face[3]), (face[3],face[1])];
		for (p1, p2) in edges
			boundary[lidx] = ( p1, p2 )
			for bidx in 1:lidx-1
				boundary_edge = boundary[bidx]
				if (p1,p2) == boundary_edge || (p2,p1) == boundary_edge
					isRepeated[lidx] = true;
					isRepeated[bidx] = true;
					nrep += 2;
				end
			end # each previous boundary --> check that the current edge has not already been added
		lidx += 1
		end # each of the 3 edges --> add it to "boundary" and check if it is not repeated
	end # each visible face --> go through its 3 edges

	"""
		2-. From each unique edge in the boundary -> build a new face with the correct neighbouring data
	"""
	numNewFaces = lidx - 1 - nrep;

	change = numNewFaces - numVisibleFaces;

	#( numNewFaces < numVisibleFaces ) && ( println( "Delete faces" ) );

	lidx = 0;
	for bidx in 1:length(boundary)
		if !(isRepeated[bidx])
			lidx += 1
			newfaceidx = ( lidx <= numVisibleFaces ) ? visibleFaces[lidx] : lastFace + lidx - numVisibleFaces;
			faces[newfaceidx] = ( boundary[bidx][1], boundary[bidx][2], newpointindex );
			norms[newfaceidx] = faceNormalCenterSIMP( points, ( boundary[bidx][1], boundary[bidx][2], newpointindex ), center )
		end
	end


	if ( change < 0 )
		cont = 0
		for start in visibleFaces[numVisibleFaces+change+1:numVisibleFaces]
			for fidx in start-cont:lastFace-cont
				faces[fidx] = faces[fidx+1]
				norms[fidx] = norms[fidx+1]
			end
			cont += 1;
		end
	end

	return lastFace + change
end
