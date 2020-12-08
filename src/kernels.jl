# short named functions, not flexible, bad assumptions or improper namings might be commited

function crossCorrelate( image::Array{T1,2}, kernel::Array{T2,2}; typ=Float32 ) where {T1<:Real,T2<:Real}
	return paddedCrossCorrelation( image, kernel, typ=typ ); 
end

function smooth( image::Array{T1,2}, size ) where {T1<:Real}
	meanKern = meanKernel( size ); 
	return crossCorrelate( image, meanKern )
end

function meanKernel( size::Integer; typ=Float32 )
	return ones( typ, size, size )./( size*size )
end

# separable edge kernels
function Prewitt_h( ; typ=Float32 )
	p1, z0, n1 = typ( 1 ), typ( 0 ), typ( -1 )
	return [ p1 z0 n1; p1 z0 n1; p1 z0 n1 ]
end

function Prewitt_v( ; typ=Float32 )
	p1, z0, n1 = typ( 1 ), typ( 0 ), typ( -1 )
	return [ p1 p1 p1; z0 z0 z0; n1 n1 n1 ]
end

function paddedCrossCorrelation( image::Array{T1,2}, kernel::Array{T2,2}; typ=Float32 ) where {T1<:Real,T2<:Real}

	N1, N2 = size( image  ); 
	M1, M2 = size( kernel ); 

	paddedImage = similar(image, N1+2*M1-2, N2+2*M2-2); 
	paddedImage[ M1:(M1+N1-1), M2:(M2+N2-1) ] .= image; 

	corrmatrix = zeros( typ, N1 + M1 - 1, N2 + M2 - 1 )

	for col in 1:( N2 + M2 - 1)
		for row in 1:( N1 + M1 - 1 )
			corrmatrix[ row, col ] = sum( paddedImage[ row:row+M1-1, col:col+M2-1 ] .* kernel )
		end
	end

	return corrmatrix
end