
function convexHull( points::Array{Tuple{T,T,T},1} ) where {T<:Real}

    @assert length(points) >= 4
    # length(points) for now, to have a static size

    # Each face consistes of 3 Vertices and 3 edges. 
	# The edges of a face, face[idx], are given by face[1]-face[2], face[2]-face[3], face[3]-face[1];
    faces = Array{Tuple{UInt16,UInt16,UInt16},1}(undef,length(points)*2);

    # Each face is in contact with 3 more faces
    neigh = Array{Tuple{UInt16,UInt16,UInt16},1}(undef,length(points)*2);

    # One normal for each face
    norms = Array{Tuple{Float32,Float32,Float32},1}(undef,length(points)*2);

    # Intializing first 4 faces and neighbouring info NOTE: need to assert that 4 first points are not colinear/coplanar
    faces[1] = ( 1, 2, 3 ); neigh[1] = ( 4, 2, 3 ); # edge 1-2 is shared with faces[4], edge 2-3 with faces[2], and edge 3-1 with faces[3]
    faces[2] = ( 2, 3, 4 ); neigh[2] = ( 1, 3, 4 ); # edge 2-3 is shared with faces[1], edge 3-4 with faces[3], and edge 4-2 with faces[4]
    faces[3] = ( 3, 4, 1 ); neigh[3] = ( 2, 4, 1 ); # edge 3-4 is shared with faces[2], edge 4-1 with faces[4], and edge 1-3 with faces[1]
    faces[4] = ( 4, 1, 2 ); neigh[4] = ( 3, 1, 2 ); # edge 4-1 is shared with faces[3], edge 1-2 with faces[1], and edge 2-4 with faces[2]

    # Faces will be modified and added. When adding new faces, "last" keeps track of the index of the last face.
    lastFace = 4;

	center = ( points[1] .+ points[2] .+ points[3] .+ points[4] ) ./ 4;

    areCoplanar( points[1], points[2], points[3], points[4] )  && (println("initial are coplanar"))

    for fidx in 1:4
		norms[fidx] = faceNormalCenter( points, faces[fidx], center )
    end

	# Iteratively modifying the current convex Hull

	for pidx in 5:length(points);    

		@inbounds new_point = points[pidx]; 
	
		# 1-. Checking faces visible to new_point.
		visibleFaces = Array{UInt16,1}(undef,lastFace); 
		numVisibleFaces = UInt16(0); 

		for fidx in 1:lastFace; 

			areCoplanar( points[faces[fidx][1]], points[faces[fidx][2]], points[faces[fidx][3]], new_point ) && (println("middle are coplanar")) 
 
			@inbounds vp = new_point .- points[faces[fidx][3]]; 
			@inbounds tf =  sign( sum( vp .* norms[ fidx ] ) ); 
			if tf == 1
				numVisibleFaces += UInt16(1); 
				visibleFaces[numVisibleFaces] = fidx; 
			end
		end

		lastFace = addNewFaces!( points, faces, neigh, norms, center, lastFace, 
						         pidx, visibleFaces, numVisibleFaces ); 

	end # each input point

    return faces[1:lastFace], neigh[1:lastFace], norms[1:lastFace]
end


"""
	needs to be tested --> tested in TEST 1
"""
@inline function faceNormalCenter( points, face, ctr )
    v =  points[face[2]] .- points[face[1]];
    w =  points[face[3]] .- points[face[1]];
	n = ( v[2]*w[3] - v[3]*w[2], v[3]*w[1] - v[1]*w[3], v[1]*w[2] - v[2]*w[1] )
	mag = sqrt( sum( n .* n ) ); 
	n = n ./ mag
	p = ctr .- points[face[1]]; 
	return ( sign( sum( p .* n ) ) == -1 ) ? n : -1 .* n; 
end


"""
	needs to be tested
"""
@inline function addNewFaces!( points,
                               faces::Array{Tuple{UInt16,UInt16,UInt16},1},
                               neigh::Array{Tuple{UInt16,UInt16,UInt16},1}, 
							   norms::Array{Tuple{Float32,Float32,Float32},1},
							   center, 
							   lastFace, 
							   newpointindex,
							   visibleFaces::Array{UInt16,1}, 
                               numVisibleFaces
                             ) 
	""" 
		1-. Collecting all edges from VisibleFaces in "boundary" and labelling those that are Repeated (not unique)
	"""
	boundary   = Array{Tuple{UInt16,UInt16},1}(undef, numVisibleFaces*3);
	isRepeated = zeros( Bool, numVisibleFaces*3 );

	lidx = 1; # linear index, it goes from 1 to numVisibleFaces*3
	nrep = 0; # number of repeats

	for face in faces[ visibleFaces[1:numVisibleFaces] ]
		@inbounds edges = [(face[1],face[2]), (face[2],face[3]), (face[3],face[1])];
		for (p1, p2) in edges
			@inbounds boundary[lidx] = ( p1, p2 )
			for bidx in 1:lidx-1
				@inbounds boundary_edge = boundary[bidx]
				if (p1,p2) == boundary_edge || (p2,p1) == boundary_edge # && !isRepeated[bidx] # each edge can only be shared by two faces, so this should not be necessary
					isRepeated[lidx] = true;
					isRepeated[bidx] = true; 
				end
			end # each previous boundary --> check that the current edge has not already been added
		lidx += 1
		end # each of the 3 edges --> add it to "boundary" and check if it is not repeated
	end # each visible face --> go through its 3 edges

	"""
		2-. From each unique edge in the boundary -> build a new face with the correct neighbouring data
	"""
	# NOTE: lidx is added 1 too much, thus "lidx - 1" below
	nrep = sum( isRepeated ); # nrep += 2 under isRepeated[bidx]
	numNewFaces = lidx - 1 - nrep;

	( numNewFaces < numVisibleFaces ) && ( println( "I Predict shit " ) ); 

	newfaces = Array{Tuple{UInt16,UInt16,UInt16},1}(undef,numNewFaces); 
	newneigh = Array{Tuple{UInt16,UInt16,UInt16},1}(undef,numNewFaces);

	lidx = 1; # linear index, from 1 to numVisibleFaces*3
	nidx = 1; # new face index, from 1 to numNewFaces

	for fidx in visibleFaces[1:numVisibleFaces]
		@inbounds face  = faces[fidx];
		@inbounds edges = [(face[1],face[2]), (face[2],face[3]), (face[3],face[1])];
		for eidx in 1:3

			@inbounds edge = edges[eidx]; 

			if !isRepeated[lidx]
				# here
				newfaces[nidx] = ( edge[1], edge[2], newpointindex );
				newneigh[nidx] = ( neigh[fidx][eidx], UInt16(0), UInt16(0) );

				newfaceposition = ( nidx <= numVisibleFaces ) ? visibleFaces[nidx] : lastFace + nidx - numVisibleFaces;

				for oidx in 1:nidx-1

					otherface = newfaces[oidx];
					otherfaceposition = ( oidx <= numVisibleFaces ) ? visibleFaces[oidx] : lastFace + oidx - numVisibleFaces;

					# I am sure this can be done smarter
					any( edge .== otherface[1] ) && ( newneigh[ oidx ] = ( newneigh[oidx][1], newneigh[oidx][2],  newfaceposition ) )
					any( edge .== otherface[2] ) && ( newneigh[ oidx ] = ( newneigh[oidx][1],  newfaceposition, newneigh[oidx][3] ) )
					
				    any( otherface .== edge[1] ) && ( newneigh[nidx] = ( newneigh[nidx][1], newneigh[nidx][2], otherfaceposition ) )
					any( otherface .== edge[2] ) && ( newneigh[nidx] = ( newneigh[nidx][1], otherfaceposition, newneigh[nidx][3] ) )

				end # for other new faces --> add neighbouring data

				nidx += 1;

			end	# for each unique edge --> build a face between that edge and the newpoint
		
		lidx += 1
	
		end # for each edge --> check that it is unique
	end # for each visible Face --> go over its 3 edges


	"""
		Either replacing or pushing the new faces with their neighbouring data to "faces" and "neigh"
	"""
	for nidx in 1:numNewFaces
		newfaceposition = ( nidx <= numVisibleFaces ) ? visibleFaces[nidx] : lastFace + nidx - numVisibleFaces;
		# here
		faces[ newfaceposition ] = newfaces[ nidx ]
		neigh[ newfaceposition ] = newneigh[ nidx ]
		norms[ newfaceposition ] = faceNormalCenter( points, newfaces[ nidx ], center )
	end


	if ( numNewFaces < numVisibleFaces ) 
		tobeDeleted = numVisibleFaces - numNewFaces
		for start in visibleFaces[1+numVisibleFaces-tobeDeleted:numVisibleFaces]
			for fidx in start:lastFace
				faces[fidx] = faces[fidx+1]
			end			
		end

		lastFace -= tobeDeleted
		# I still need to fix neighbourhoods
	end

	# C
	#println( lastFace, "  ", numVisibleFaces, "  ", nrep, "  ", numNewFaces, "  ", lastFace + numNewFaces - numVisibleFaces )
	# return new lastFace
	change = ( numNewFaces < numVisibleFaces ) ? 0 : numNewFaces - numVisibleFaces; 
	return lastFace + change
end

function areCoplanar( x, y, z, p )

    a1 = x[2] - x[1] ; 
    b1 = y[2] - y[1] ; 
    c1 = z[2] - z[1] ; 
    a2 = x[3] - x[1] ; 
    b2 = y[3] - y[1] ; 
    c2 = z[3] - z[1] ; 
    a = b1 * c2 - b2 * c1 ; 
    b = a2 * c1 - a1 * c2 ; 
    c = a1 * b2 - b1 * a2 ; 
    d = - a * x[1] - b * y[1] - c * z[1];     

    return a * p[1] + b * p[2] + c * p[3] + d == 0 
end
