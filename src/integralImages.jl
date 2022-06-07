"""
	2D
"""

function integralArray( img::AbstractArray{<:Real,2}; type=Float32, fun=(x)->(x) )
	return integralArray!( img, zeros( type, size(img) .+ 1 ), fun=fun )
end

function integralArray!( img::AbstractArray{<:Real,2}, intArr::AbstractArray{<:AbstractFloat,2}; fun=(x)->(x) ) 
	
	@inbounds for c in 1:size(img,2), r in 1:size(img,1)
		intArr[1+r,1+c] = fun(img[r,c]) + intArr[1+r-1,1+c] + intArr[1+r,1+c-1] - intArr[1+r-1,1+c-1]
	end
	return intArr
end

function integralArea( intArr::AbstractArray{<:Real,2}, TL, BR )
	TL   = TL .+ 1;
	BR   = BR .+ 1;
	area = intArr[BR[1],BR[2]] - intArr[BR[1],TL[2]] - intArr[TL[1],BR[2]] + intArr[TL[1],TL[2]]
end

# integral array for ZNCC, copied from quickPIV

function integralArrayZNCC( padData::AbstractArray{Complex{T},N} ) where {T<:AbstractFloat,N}
	return integralArrayZNCC!( padData, zeros( T, size(padData).+1 ) ); 
end

function integralArrayZNCC!( padData::AbstractArray{Complex{T},2}, intArr::Array{T,2} ) where {T<:AbstractFloat}

	h, w = size(intArr) .- 1

	@inbounds for c in 2:w+1, r in 2:h+1
		intArr[r,c] = (padData[r-1,c-1].re)^2 + intArr[r-1,c] + intArr[r,c-1] - intArr[r-1,c-1]
	end

	return intArr
end

function integralArrayZNCC!( padData::AbstractArray{Complex{T},3}, intArr::Array{T,3} ) where {T<:AbstractFloat} 

	h, w, d = size(intArr) .- 1;
	@inbounds for z in 2:d+1, c in 2:w+1
		tmp = 0.0; 
		for r in 2:h+1      
			val2          = real( padData[r-1,c-1,z-1] )^2
			intArr[r,c,z] = val2 + intArr[r,c-1,z] + intArr[r,c,z-1] - intArr[r,c-1,z-1] + tmp;
			tmp          += val2; 
	end end

	return intArr
end

# what did IIM mean...?

function integralArrayIIM!( img::AbstractArray{<:Real,2}, intArr::AbstractArray{<:AbstractFloat,3}, z ) 
	
	@inbounds for c in 1:size(img,2), r in 1:size(img,1)
		intArr[1+r,1+c,z] = img[r,c] + intArr[1+r-1,1+c,z] + intArr[1+r,1+c-1,z] - intArr[1+r-1,1+c-1,z]
	end
	return intArr
end

function integralAreaIIM( intArr::AbstractArray{<:Real,3}, z,TL, BR )
	TL   = TL .+ 1;
	BR   = BR .+ 1;
	area = intArr[BR[1],BR[2],z] - intArr[BR[1],TL[2],z] - intArr[TL[1],BR[2],z] + intArr[TL[1],TL[2],z]
end

# Considering offset of the img within the integral array

function integralArray_subimg( img::AbstractArray{<:Real,2}, range=(0,0,0,0); type=Float32, fun=(x)->(x) )
	return integralArray_subimg!( img, zeros( type, size(img) .+ 1 ), off, fun=fun )
end

function integralArray_subimg!( img::AbstractArray{<:Real,2}, intArr::AbstractArray{<:AbstractFloat,2}, range=(0,0,0,0); fun=(x)->(x)) 
	
	ymin, ymax = range[1], range[2]; 
	xmin, xmax = range[3], range[4];

	@inbounds for c in xmin:xmax, r in ymin:ymax
		ir, ic = r - ymin, c - xmin; 
		intArr[1+ir,1+ic] = fun(img[r,c]) + intArr[1+ir-1+off[1],1+ic] + intArr[1+ir,1+ic-1] - intArr[1+ir-1,1+ic-1]
	end
	return intArr
end

# Considering offset of the img within the integral array

function integralArray_off( img::AbstractArray{<:Real,2}, off=(0,0); type=Float32, fun=(x)->(x) )
	return integralArray_off!( img, zeros( type, size(img) .+ 1 ), off, fun=fun )
end

function integralArray_off!( img::AbstractArray{<:Real,2}, intArr::AbstractArray{<:AbstractFloat,2}, off=(0,0); fun=(x)->(x)) 
	
	@inbounds for c in 1:size(img,2), r in 1:size(img,1)
		intArr[1+r+off[1],1+c+off[2]] = fun(img[r,c]) + intArr[1+r-1+off[1],1+c+off[2]] + intArr[1+r+off[1],1+c-1+off[2]] - intArr[1+r-1+off[1],1+c-1+off[2]]
	end
	return intArr
end

# Specific function to compute mean maxima and minima, which requires padded integral arrays

function integralArray_pad( img::AbstractArray{<:Real,2}, pad=(0,0); type=Float32, fun=(x)->(x) )
	return integralArray_pad!( img, zeros( type, size(img) .+ 1 .+ 2 .* pad ), pad, fun=fun )
end

function integralArray_pad!( img::AbstractArray{<:Real,2}, intArr::AbstractArray{<:AbstractFloat,2}, pad; fun=(x)->(x) ) 
	
	for c in 1:size(img,2)      # c    = column in input image
		c_ia = 1+c+pad[2];      # c_ia = offset column index in integral array
		for r in 1:size(img,1)  # r    = row in input image
			r_ia = 1+r+pad[1];  # r_ia = offset row index in integral array
			intArr[r_ia,c_ia] = fun(img[r,c]) + intArr[r_ia-1,c_ia] + intArr[r_ia,c_ia-1]- intArr[r_ia-1,c_ia-1]
	end end

	lastr = 1+size(img,1)+pad[1]; 
	lastc = 1+size(img,2)+pad[2];

	for c in 1:lastc
		val = intArr[lastr,c]
		for r in lastr+1:size(intArr,1)
			intArr[r,c] = val
	end end

	for c in lastc+1:size(intArr,2)
		intArr[:,c] .= intArr[:,c-1]
	end

	return intArr
end

function integralArray_pad( vol::AbstractArray{<:Real,3}, pad=(0,0,0); type=Float32, fun=(x)->(x) )
	return integralArray_pad!( vol, zeros( type, size(vol) .+ 1 .+ 2 .*pad ), pad, fun=fun )
end

function integralArray_pad!( vol::AbstractArray{<:Real,3}, intArr::AbstractArray{<:AbstractFloat,3}, pad; fun=(x)->(x)) 

	for z in 1:size(vol,3)         # z    = depth index in input volume
		z_ia = z+pad[3]+1;         # z_ia = offset depth index in integral array
		for c in 1:size(vol,2)     # c    = columns in input volume
			c_ia = c+pad[2]+1;     # c_ia = offset column index in integral array
			reuse = 0.0
			for r in 1:size(vol,1) # r    = row in input volume
				r_ia = r+pad[1]+1; # r_ia = offset row index in integral array
				tmp = fun(vol[r,c,z])
				intArr[r_ia,c_ia,z_ia] = tmp + intArr[r_ia,c_ia-1,z_ia] + intArr[r_ia,c_ia,z_ia-1] - intArr[r_ia,c_ia-1,z_ia-1] + reuse; 
				reuse += tmp
	end end end

	lastr = 1+size(vol,1)+pad[1]; 
	lastc = 1+size(vol,2)+pad[2];
	lastz = 1+size(vol,3)+pad[3];

	for z in 1:lastz
		for c in 1:lastc
			val = intArr[lastr,c,z]
			for r in lastr+1:size(intArr,1)
				intArr[r,c,z] = val
		end end 

		for c in lastc+1:size(intArr,2)
			intArr[:,c,z] .= intArr[:,c-1,z]
		end
	end

	for z in lastz+1:size(intArr,3)
		intArr[:,:,z] .= intArr[:,:,z-1]
	end

	return intArr
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


function integralArray( vol::AbstractArray{<:Real,3}; type=Float32, fun=(x)->(x) )

	return integralArray!( vol, zeros( type, size(vol) .+ 1 ), fun )
end

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

function integralArray_masked( vol::AbstractArray{<:Real,3}, mask; type=Float32, fun=(x)->(x) )

	return integralArray_masked!( vol, mask, zeros( type, size(vol) .+ 1 ), fun )
end


function integralArray_masked!( vol::AbstractArray{<:Real,3}, mask, intArr::AbstractArray{<:AbstractFloat,3}, fun=(x)->(x) )

	@inbounds for z in 2:size(vol,3)+1, c in 2:size(vol,2)+1
		reuse = 0.0
		for r in 2:size(vol,1)+1
			tmp = fun(vol[r-1,c-1,z-1])*mask[r-1,c-1,z-1]
			intArr[r,c,z] = tmp + intArr[r,c-1,z] + intArr[r,c,z-1] - intArr[r,c-1,z-1] + reuse; 
			reuse += tmp
		end
	end
	return intArr
end

function integralArray_subvol( img::AbstractArray{<:Real,3}, range=(0,0,0,0,0,0); type=Float32, fun=(x)->(x) )
	return integralArray_subvol!( img, zeros( type, size(img) .+ 1 ), range, fun=fun )
end

function integralArray_subvol!( img::AbstractArray{<:Real,3}, intArr::AbstractArray{<:AbstractFloat,3}, range=(0,0,0,0,0,0); fun=(x)->(x)) 
	
	ymin, ymax = range[1], range[2]; 
	xmin, xmax = range[3], range[4];
	zmin, zmax = range[5], range[6];

	@inbounds for z in zmin:zmax, c in xmin:xmax
		iz, ic = z - zmin + 1, c - xmin + 1; 
		reuse = 0.0
		for r in ymin:ymax
			ir = r - ymin + 1;
			tmp = fun(img[r,c,z])
			intArr[1+ir,1+ic,1+iz] = tmp + intArr[1+ir,1+ic-1,1+iz] + intArr[1+ir,1+ic,1+iz-1] - intArr[1+ir,1+ic-1,1+iz-1] + reuse; 
			reuse += tmp
		end
	end
	return intArr
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

