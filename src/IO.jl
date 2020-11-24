function grayscale255( filename::String; path::String=pwd() ) 

	img = Float64.( ColorTypes.Gray.( load( path*filename ) ) )
	mn, mx = extrema( img )

	im255 = Array{ Int64, 2 }(undef, size(img) )
	for e in 1:length( img ) 
		im255[e] = 1 + round( Int64, ( img[e] - mn )/( mx - mn ) * 254 )
	end

	return im255
end

	
	
