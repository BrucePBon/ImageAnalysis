function mean( img::Array{T,2} ) where {T<:Real}

	ret = 0.0
	N = length( img ); 

	@simd for li in 1:N
		@inbounds ret += img[li]
	end

	return ret/N
end

function variance( img::Array{T,2} ) where {T<:Real}

	mean = mean( img )

	ret = 0.0
	N = length( img ) 

	@simd for li in 1:N
		@inbounds ret += ( img[li] - mean )*( img[li] - mean )
	end

	return ret/N
end

function variance( img::Array{T,2}, mean ) where {T<:Real}

	ret = 0.0
	N = length( img ) 

	@simd for li in 1:N
		@inbounds ret += ( img[li] - mean )*( img[li] - mean )
	end

	return ret/N
end

function std( img::Array{T,2} ) where {T<:Real}

	return sqrt( variance( img ) )
end

function std( img::Array{T,2}, mean ) where {T<:Real}

	return sqrt( variance( img, mean ) )
end

function cov( img1::Array{T1,2}, img2::Array{T2,2}) where {T1<:Real,T2<:Real}

	length( img1 ) !== length( img2 ) && ( error("Arrays have different lengths") )
	
	m1 = mean( img1 ) 
	s1 = std( img1, m1 ) 

	m2 = mean( img2 )
	s2 = std( img2, m2 ) 

	ret = 0.0
	N = length( img1 ) 

	@simd for li in 1:N 
		@inbounds ret += ( img1[li] - m1 )*( img2[li] - m2 ) 
	end

	return ret/N
end

function correlation( img1::Array{T1,2}, img2::Array{T2,2} ) where {T1<:Real,T2<:Real}

	length( img1 ) !== length( img2 ) && ( error("Arrays have different lengths") )
	
	m1 = mean( img1 ) 
	s1 = std( img1, m1 ) 

	m2 = mean( img2 )
	s2 = std( img2, m2 ) 

	ret = 0.0
	N = length( img1 ) 

	@simd for li in 1:N 
		@inbounds ret += ( img1[li] - m1 )*( img2[li] - m2 ) 
	end

	return ret/(N*s1*s2)
end