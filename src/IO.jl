function grayscale( filename::String; path::String=pwd(), ftyp=Float32, ityp=Int16 )

	occursin( '/', filename ) && ( path = "" )

	img = FileIO.load( path*filename )

	# Temporary array, used to store the grayscale values of "img" as Floats

	imgf = zeros( ftyp, size( img ) )

	for idx in 1:length(img)
		imgf[idx] = ftyp( ColorTypes.Gray( img[idx] ) )
	end

	# The grayscale image is stored in "imgG" as Integer values

	imgG = Array{ ityp, 2 }(undef, size(img) )
	for li in 1:length( imgG )
		imgG[li] = 1 + round( ityp, imgf[li]* 254 )
	end

	return imgG
end

function grayscale!( filename::String, img::Array{<:Integer,2}; path::String=pwd(), ftyp=Float32 )

	occursin( '/', filename ) && ( path = "" )
	( isfile( path*filename) ) || error( "file does not exists" )

	img_ = FileIO.load( path*filename )

	( size( img_ ) == size(img) ) || error( "input image should have dimensions: ", size(img_) )

	T = eltype( img );

	for li in 1:length( img )
		img[li] = 1 + round( T, ftyp( ( img_[li].val ) )*254 )
	end
end


function grayscale255( filename::String; path::String=pwd(), typ=Float32 )

	occursin( '/', filename ) && ( path = "" )

	img = FileIO.load( path*filename )

	imgf = zeros( typ, size( img ) )

	for idx in 1:length(img)
		imgf[idx] = typ( ColorTypes.Gray( img[idx] ) )
	end

	mn, mx = extrema( imgf )

	im255 = Array{ Int64, 2 }(undef, size(img) )
	for e in 1:length( img )
		im255[e] = 1 + round( Int64, ( imgf[e] - mn )/( mx - mn ) * 254 )
	end

	return im255
end

function grayscale255rev( filename::String; path::String=pwd() )

	occursin( '/', filename ) && ( path = "" )

	img = load( path*filename )

	imgf = zeros( Float64, size( img ) )

	for idx in 1:length(img)
		imgf[idx] = Float64( ColorTypes.Gray( img[idx] ) )
	end

	mn, mx = extrema( imgf )

	im255 = Array{ Int64, 2 }(undef, size(imgf) )
	for e in 1:length( img )
		im255[e] = 256 - 1 - round( Int64, ( imgf[e] - mn )/( mx - mn ) * 254 )
	end

	return im255
end


function RGBfloat( filename::String; path::String=pwd(), typ=Float32 )

	occursin( '/', filename ) && ( path = "" )

	data = load( path*filename )

	RGBimg = Array{ typ, 3 }( undef, size( data, 1 ), size( data, 2 ), 3 )

	Goff = length( data );
	Boff = Goff*2;

	@inbounds for li in 1:length(data)
	    pixel = data[li]
	    RGBimg[   li  ] = typ( data[li].r )
	    RGBimg[li+Goff] = typ( data[li].g )
	    RGBimg[li+Boff] = typ( data[li].b )
	end

	return RGBimg
end
