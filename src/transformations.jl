#=
BILINEAR INTERPOLATION
	Q1      Q2

	

	Q3      Q4
=#

function bilinear_distances( size )

	distances = Array{Tuple{Float32,Float32},1}(undef,size)
	delta     = 1/(size+1); 

	d = delta
	for idx in 1:size
		distances[idx] = ( 1-d, d )
		d += delta
	end
	
	return distances
end

function discrete_bilinear( Q12, Q34, distances, id1, id2 )
	
	intp1 = sum( distances[id1] .* Q12 )
	intp2 = sum( distances[id1] .* Q34 )
	return  sum( distances[id2] .* (intp1,intp2) )
end

function bilinear_interpolation( img, grid; sz=4 )

	 h,  w = size(img); 
	gh, gw = size(grid); 

	startx, stepx = w/(gw+2)
	starty, stepy = h/(gh+2)

	for
end
