# short named functions, not flexible, bad assumptions or improper namings might be commited

function convolution( image::Array{T1,N}, kernel::Array{T2,N}; typ=Float32 ) where {T1<:Real,T2<:Real,N}
	return FFTConvolution( image, kernel, typ=typ );
end

function convolution_b( image::Array{T1,2}, kernel::Array{T2,2}; typ=Float32 ) where {T1<:Real,T2<:Real}
	return FFTConvolution_b( image, kernel, typ=typ );
end

function smooth( image::Array{T1,2}, size ) where {T1<:Real}
	if size == 0
		return image
	end
	meanKern = meanKernel( size );
	return FFTConvolution( image, meanKern )
end

function gradient( image::Array{T1,2}; kernel="Prewitt" ) where {T1<:Real}

	kh = ( kernel=="Sobel" ) ? Sobel_h() : Prewitt_h();
    kv = ( kernel=="Sobel" ) ? Sobel_v() : Prewitt_v();

	Gy = convolution( image, kh )
	Gx = convolution( image, kv )

	return Gy[2:end-1,2:end-1], Gx[2:end-1,2:end-1]
end

# smoothing kernel

function meanKernel( size::Integer; typ=Float32 )
	return ones( typ, size, size )./( size*size )
end

# edge kernels

# Prewitt, separable
function Prewitt_h( ; typ=Float32 )
	p1, z0, n1 = typ.( ( 1, 0, -1 ) )
	return [ p1 z0 n1; p1 z0 n1; p1 z0 n1 ]
end

function Prewitt_v( ; typ=Float32 )
	p1, z0, n1 = typ.( ( 1, 0, -1 ) )
	return [ p1 p1 p1; z0 z0 z0; n1 n1 n1 ]
end

# Sobel, separable
function Sobel_h( ; typ=Float32 )
	p1, p2, z0, n1, n2 = typ.( ( 1, 2, 0, -1, -2 ) )
	return [ p1 z0 n1; p2 z0 n2; p1 z0 n1 ]
end

function Sobel_v( ; typ=Float32 )
	p1, p2, z0, n1, n2 = typ.( ( 1, 2, 0, -1, -2 ) )
	return [ p1 p2 p1; z0 z0 z0; n1 n2 n1 ]
end
