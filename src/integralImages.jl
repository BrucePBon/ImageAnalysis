"""
	2D
"""

function integralArray!( img::AbstractArray{<:Real,2}, intArr::AbstractArray{<:AbstractFloat,2}, fun=(x)->(x) ) 
	
	@inbounds for c in 2:size(img,2)+1, r in 2:size(img,1)+1
		intArr[r,c] = fun(img[r-1,c-1]) + intArr[r-1,c] + intArr[r,c-1] - intArr[r-1,c-1]
	end
	return intArr
end

function integralArray( img::AbstractArray{<:Real,2}; type=Float32, fun=(x)->(x) )

	return integralArray!( img, zeros( type, size(img) .+ 1 ), fun )
end



function integralArea( intArr::AbstractArray{<:Real,2}, TL, BR )
	TL   = TL .+ 1;
	BR   = BR .+ 1;
	area = intArr[BR[1],BR[2]] - intArr[BR[1],TL[2]] - intArr[TL[1],BR[2]] + intArr[TL[1],TL[2]]
end



# Reversed variant, which is used during normalized cross-correlation
function integralArrayRev!( img::AbstractArray{<:Real,2}, intArr::AbstractArray{<:AbstractFloat,2}, fun=(x)->(x) ) 
	h, w = size(img)
	@inbounds for c in 2:w+1, r in 2:h+1
		intArr[r,c] = fun(img[h-(r-2),w-(c-2)]) + intArr[r-1,c] + intArr[r,c-1] - intArr[r-1,c-1]
	end
	return intArr
end

function integralArrayRev( img::AbstractArray{<:Real,2}; type=Float32, fun=(x)->(x) )

	return integralArrayRev!( img, zeros( type, size(img) .+ 1 ), fun )
end





"""
	3D: optimized integral array by reusing previous operations
"""
# Below is the unoptimized code, in case I want to check again: 
#    for z in 2:d+1, c in 2:w+1, r in 2:h+1      
#        arr[r,c,z] = fun(vol[r-1,c-1,z-1]) + arr[r-1,c,z] + arr[r,c-1,z] + arr[r,c,z-1] -
#                     arr[r,c-1,z-1] - arr[r-1,c,z-1] - arr[r-1,c-1,z] + arr[r-1,c-1,z-1]
#    end

function integralArray!( vol::AbstractArray{<:Real,3}, intArr::AbstractArray{<:AbstractFloat,3}, fun=(x)->(x) )

	@inbounds for z in 2:size(vol,3)+1, c in 2:size(vol,2)+1
		reuse = 0.0
		for r in 2:size(vol,1)+1
			tmp = fun(vol[r-1,c-1,z-1])
			intArr[r,c,z] = tmp + intArr[r,c-1,z] + intArr[r,c,z-1] - intArr[r,c-1,z-1] + reuse; 
			reuse += tmp
		end
	end
	return intArr
end

function integralArray( vol::AbstractArray{<:Real,3}; type=Float32, fun=(x)->(x) )

	return integralArray!( vol, zeros( type, size(vol) .+ 1 ), fun )
end



function integralArea( intArr::AbstractArray{<:Real,3}, TL, BR )
	TL = TL .+ 1; 
	BR = BR .+ 1; 
	area  = intArr[BR[1],BR[2],BR[3]] - intArr[TL[1],TL[2],TL[3]]
	area -= intArr[TL[1],BR[2],BR[3]] + intArr[BR[1],TL[2],BR[3]] + intArr[BR[1],BR[2],TL[3]]
	area += intArr[BR[1],TL[2],TL[3]] + intArr[TL[1],BR[2],TL[3]] + intArr[TL[1],TL[2],BR[3]]
    return area
end



# Reversed variant, which is used during normalized cross-correlation
function integralArrayRev!( vol::AbstractArray{<:Real,3}, intArr::AbstractArray{<:AbstractFloat,3}, fun=(x)->(x) )

	h, w, d = size(vol)
	@inbounds for z in 2:d+1, c in 2:w+1
		reuse = 0.0
		for r in 2:h+1
			tmp = fun(vol[h-(r-2),w-(c-2),d-(z-2)])
			intArr[r,c,z] = tmp + intArr[r,c-1,z] + intArr[r,c,z-1] - intArr[r,c-1,z-1] + reuse; 
			reuse += tmp
		end
	end
	return intArr
end

function integralArrayRev( vol::AbstractArray{<:Real,3}; type=Float32, fun=(x)->(x) )

	return integralArrayRev!( vol, zeros( type, size(vol) .+ 1 ), fun )
end





"""
	Derived functions
"""
function scan4Simmetry( intArr, rad ) 
	
	score = zeros( length( rad+1:size(intArr,1)-rad), length( rad+1:size(intArr,2)-rad ) )

	cont = 1; 
	for c in rad+1:size(intArr,2)-rad;     c1, c2 = c - rad, c + rad
		for r in rad+1:size(intArr,1)-rad; r1, r2 = r - rad, r + rad
		
			Q1  = integralArea( intArr, (r1,c ), ( r-1,c2-1) )
			Q2  = integralArea( intArr, (r1,c1), ( r-1,c -1) )
			Q3  = integralArea( intArr, (r ,c1), (r2-1,c -1) )
			Q4  = integralArea( intArr, (r ,c ), (r2-1,c2-1) )

			acc = Q1 + Q2 + Q3 + Q4; 
			ps  = [ Q1, Q2, Q3, Q4 ] ./ acc 

			score[ cont ] = 1/( 1 + sum( abs.( ps .- 0.25 ) ) )
			score[ cont ] = isnan( score[cont] ) ? 0.0 : score[cont]; 
			cont += 1
	end	end

	return score	
end

function followCentroids( data, intArr, rad; iters=5 ) 

	# Goal: At any position, compute centroid, move to centroid
	# recompute centroid and so on until converging

	score = zeros( Int64, size(data) )

	xx = [ x for y in -rad:rad, x in -rad:rad ]
	yy = [ y for y in -rad:rad, x in -rad:rad ]

	h, w = size( data ); 
	for c in 1:w;    
		for r in 1:h; 

			pos = r, c; 
			#println( pos )
			# Contemplating the case where the position is closer than "rad"
			# to the edges of the images. As we iterate and update the position
			# of the centroid, it becomes necessary. Here is not needed though.
			left   = rad + min(0, pos[1]-1-rad)
			up     = rad + min(0, pos[2]-1-rad)
			right  = min(  h - pos[ 1 ], rad  )
			down   = min(  w - pos[ 2 ], rad  )
		
			TL     = pos .- ( left , up   ) .- 1
			BR     = pos .+ ( right, down )
		
			acc    = integralArea( intArr, TL, BR )
			( acc == 0 ) && ( continue; )

			subimg = pos[1]-left:pos[1]+right, pos[2]-up:pos[2]+down
			subxy  =  rad+1-left:rad+1+right ,  rad+1-up:rad+1+down
		
			xmean  = sum( (data[ subimg... ] ./ acc) .* xx[ subxy... ] )
			ymean  = sum( (data[ subimg... ] ./ acc) .* yy[ subxy... ] )

			for i in 1:iters

				offset = round.( Int64, (ymean,xmean) )

				if offset[1] == 0 && offset[2] == 0
					break
				end

				pos    = pos .+ offset
				#println( "\t", pos )
				left   = rad + min(0, pos[1]-1-rad)
				up     = rad + min(0, pos[2]-1-rad)
				right  = min(  h - pos[ 1 ], rad  )
				down   = min(  w - pos[ 2 ], rad  )
			
				TL     = pos .- ( left , up   ) .- 1
				BR     = pos .+ ( right, down )
			
				acc    = integralArea( intArr, TL, BR )
				( acc == 0 ) && ( break; )

				subimg = pos[1]-left:pos[1]+right, pos[2]-up:pos[2]+down
				subxy  =  rad+1-left:rad+1+right ,  rad+1-up:rad+1+down
			
				xmean  = sum( (data[ subimg... ] ./ acc) .* xx[ subxy... ] )
				ymean  = sum( (data[ subimg... ] ./ acc) .* yy[ subxy... ] )
			end

			pos = pos .+ round.( Int64, (ymean,xmean) )

			score[ pos... ] += 1; 
	end	end

	return score
end

function integralAreaZNCC( sumS::Array{T,2}, sumS2::Array{T,2},
                           row, col, ss, si ) where {T<:Real}
    pad = 1;
    row = row + pad;
    col = col + pad;
    r0, r1 = max( 1, row - si[1] ), min( ss[1] + pad, row );
    c0, c1 = max( 1, col - si[2] ), min( ss[2] + pad, col );

    opS  = 0.0
    opS += sumS[ r1, c1 ];
    opS -= sumS[ r1, c0 ];
    opS -= sumS[ r0, c1 ];
    opS += sumS[ r0, c0 ];

    opS2  = 0.0
    opS2 += sumS2[ r1, c1 ];
    opS2 -= sumS2[ r1, c0 ];
    opS2 -= sumS2[ r0, c1 ];
    opS2 += sumS2[ r0, c0 ];

    return opS, opS2
end

