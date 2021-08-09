function convexHullSIMP!( points::Array{Tuple{T,T,T},1}; typ=UInt16 ) where {T<:Real}

    @assert length(points) >= 4
    # length(points) for now, to have a static size

    faces = Array{Tuple{typ,typ,typ},1}(undef,length(points)*2);
    norms = Array{Tuple{Float32,Float32,Float32},1}(undef,length(points)*2);
	visibleFaces = Array{typ,1}(undef,length(points)*2); 

	# making sure points are not the same, colinear or coplanar
	idx = 2
	while areSame( points[1], points[idx] ) && idx < length(points)
		idx += 1; 
	end	
	if ( idx !== 2 )
		temp = points[idx]
		points[idx] = points[2]
		points[2] = temp
	end

	idx = 3
	while areColinear( points[1], points[2], points[idx] ) && idx < length(points)
		idx += 1;
	end
	if ( idx !== 3 )
		temp = points[idx]
		points[idx] = points[3]
		points[3] = temp
	end

	idx = 4
	while areCoplanar( points[1], points[2], points[3], points[idx] ) && idx < length(points)
		idx += 1;
	end
	if ( idx !== 4 )
		temp = points[idx]
		points[idx] = points[4]
		points[4] = temp
	end
	
    faces[1] = ( 1, 2, 3 );
    faces[2] = ( 2, 3, 4 );
    faces[3] = ( 3, 4, 1 );
    faces[4] = ( 4, 1, 2 );

    lastFace = 4;
	center   = ( points[1] .+ points[2] .+ points[3] .+ points[4] ) ./ 4;
    for fidx in 1:4
		norms[fidx] = faceNormalCenterSIMP( points, faces[fidx], center )
    end

	@inbounds for pidx in 5:length(points);  
  
		new_point = points[pidx]; 

		numVisibleFaces = typ(0); 

		for fidx in 1:lastFace; 
			
			vp = new_point .- points[faces[fidx][1]]; 

			if sign( sum( vp .* norms[ fidx ] ) ) == 1
				numVisibleFaces += typ(1); 
				visibleFaces[numVisibleFaces] = fidx; 
			end
		end

		lastFace = addNewFacesSIMP!( points, faces, norms, center, lastFace, 
						             pidx, visibleFaces, numVisibleFaces ); 

	end # each input point

    return faces[1:lastFace], norms[1:lastFace]
end

"""
	needs to be tested --> tested in TEST 1
"""
# not normalize
@inline function faceNormalCenterSIMP( points, face, ctr )
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
# can this be made more compact?
@inline function addNewFacesSIMP!( points,
                               faces::Array{Tuple{T,T,T},1},
							   norms::Array{Tuple{Float32,Float32,Float32},1},
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
	change      = numNewFaces - numVisibleFaces;

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

cross( a, b ) = a[2]*b[3] - a[3]*b[2], a[3]*b[1] - a[1]*b[3], a[1]*b[2] - a[2]*b[1]

function areSame( x, y )
 
	#println( "same test: ", (x, y) )  
	return x == y
end

function areColinear( x, y, z )
	
	v1 = x .- z
	v2 = y .- z
	n  = cross( v1, v2 )

	#println( "colinear test: ", n )

	return n[1] == 0 && n[2] == 0 && n[3] == 0
end

function areCoplanar( x, y, z, p )

	v1 = ( x[1] - p[1], x[2] - p[2], x[3] - p[3] )
	v2 = ( y[1] - p[1], y[2] - p[2], y[3] - p[3] )
	v3 = ( z[1] - p[1], z[2] - p[2], z[3] - p[3] )

 	test = sum( v1 .* cross( v2, v3 ) )

	#println( "coplanar test: ", test )  

    return test == 0 
end






























function convexHullSIMP2( points::Array{Tuple{T,T,T},1}; typ=UInt16 ) where {T<:Real}

    @assert length(points) >= 4
    # length(points) for now, to have a static size

    faces = Array{Tuple{typ,typ,typ},1}(undef,0);
    norms = Array{Tuple{Float32,Float32,Float32},1}(undef,0);


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
							   norms::Array{Tuple{Float32,Float32,Float32},1},
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


function areCoplanarSIMP( x, y, z, p )

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
