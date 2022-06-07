
# Otsu thresholding
# Does it only work on integers...?

function o_varB( k, sumP, sumIP, meanT )
	w0 = sumP[k]; 
	w1 = sumP[end] - sumP[k]; 
	m0 = sumP[k]
	m1 = sumP[end] - sumP[k]; 
	return w0*w1*( m0 - m1 )^2
	return ( meanT*sumP[k] - sumIP[k] )^2/( sumP[k]*( 1 - sumP[k] ) )
end

function otsu( p::Array{Float32,1}; min=0, step=1, max=255 )

	L = length( p )

	sumP = copy( p )
	@inbounds for i in 2:L
		sumP[i] += sumP[ i - 1 ]
	end

	sumIP = zeros( Float32, L )
	sumIP[1] = min
	ith_bin  = min + step
	@inbounds for i in 2:L
		sumIP[i] = ith_bin*p[i] + sumIP[ i - 1 ]
		ith_bin += step; 
	end
	meanT = sumIP[end]

	maxK = 2
	maxvarB = o_varB( maxK, sumP, sumIP, meanT )
	for k in 3:L-1
		varB = o_varB( k, sumP, sumIP, meanT )
		if ( varB > maxvarB )
			maxK = k
			maxvarB = varB
		end
	end

	return maxK*step
end

function otsu( img, bin_data::Union{Tuple{T1,T2,T3}, AbstractRange{T4}} ) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real}

	nh, _, _ = normhistogram( img, bin_data )
	ot = otsu( nh );

	#println( "Otsu threshold = ", ot );
	mask = zeros( Bool, size(img) );

	for e in 1:length(img)
		mask[e] = img[e] > ot
	end

	return mask
end

function otsuTH( img::Array{<:Real,N}; nbins=100, min=nothing, max=nothing ) where {N}
	mn, mx = extrema( img ); 
	min = ( isnothing(min) ) ? mn : min; 
	max = ( isnothing(max) ) ? mx : max;
	bin_data = min:(max-min)/nbins:max; 
	nh, _, _ = normhistogram( img, bin_data )
	return otsu( nh )
end

function otsu( img::Array{<:Real,N} ) where {N}

	mn, mx = extrema( img ); 
	bin_data = mn:(mx-mn)/100:mx; 
	nh, _, _ = normhistogram( img, bin_data )
	ot = otsu( nh );

	#println( "Otsu threshold = ", ot );
	mask = zeros( UInt8, size(img) );

	for e in 1:length(img)
		mask[e] = reinterpret( UInt8, img[e] > ot )
	end

	return mask
end

function otsuArr( p::Array{Float32,1}; min=0, step=1, max=255 )

	arr = copy( p )
	L = length( p )

	sumP  = copy( p )
	@simd for i in 2:L
		 @inbounds sumP[i] += sumP[ i - 1 ]
	end

	sumIP = zeros( Float32, L )
	sumIP[1] = min
	ith_bin  = min + step
	@inbounds for i in 2:L
		sumIP[i] = ith_bin*p[i] + sumIP[ i - 1 ]
		ith_bin += step; 
	end
	meanT = sumIP[end]

	for k in 1:L
		arr[k] = o_varB( k, sumP, sumIP, meanT )
	end

	return arr
end
