# NOTE: the implementation in openCV is very optimal, I have to study that one and hopefully
# get to implement it.

connectivity4() = ( (-1,0),(0,1),(1,0),(0,-1) )
connectivity8() = ( (-1,0,0),(1,0,0),(0,-1,0),(0,1,0),(0,0,-1),(0,0,1) )
default_connect(N) = ( N == 2 ) ? connectivity4() : connectivity8();

""" one at a time connected component labelling in 2D """

function OAATccl( mask::Array{UInt8,2}; connectivity=connectivity4(), ref=UInt8(1) )

    h, w  = size(mask);
    mask_ = copy(mask);

    components = Array{Array{Tuple{Int64,Int64},1},1}(undef,0);
    comp_idx   =  0;
    stack      = [];

    for col in 1:w, row in 1:h

        # starting the process
        if mask_[ row, col ] == ref
            comp_idx += 1;
            push!( components, Array{Tuple{Int64,Int64},1}(undef,0) )
            push!( stack, (row,col) );
        end

        # Iterating
        while length(stack) > 0
            idx = pop!( stack )
            mask_[ idx... ] = UInt8(!Bool(ref))
            push!( components[comp_idx], idx );

            for off in connectivity
                r, c = min.((h,w), max.(1, idx .+ off))
                if mask_[r,c] == ref
                    push!( stack, ( r, c ) );
                end
            end
        end
    end
    return components;
end

function OAATcclrev( mask::Array{UInt8,2}; connectivity=connectivity4() )

	return OAATccl( mask, connectivity=connectivity, ref=UInt8(0) )
end

""" one at a time connected component labelling in OAAT 3D """

function OAATccl( mask::Array{UInt8,3}; ref=UInt8(1), connectivity=connectivity8() )

    h, w, d = size(mask);
    mask_   = copy(mask);

    components = Array{Array{Tuple{Int64,Int64,Int64},1},1}(undef,0);
    comp_idx   =  0;
    stack      = [];
	
	for zet in 1:d, col in 1:w, row in 1:h

        if mask_[ row, col, zet ] == ref
            comp_idx += 1;
            push!( components, Array{Tuple{Int64,Int64,Int64},1}(undef,0) )
            push!( stack, (row,col,zet) );
        end

        # Iterating
        while length(stack) > 0
            idx = pop!( stack )
            mask_[ idx... ] = UInt8(!Bool(ref))
            push!( components[comp_idx], idx );

            for off in connectivity
                r, c, z = min.( (h,w,d), max.( 1, idx .+ off ) )
                if mask_[r,c,z] == ref
                    push!( stack, ( r, c, z ) );
                end
            end
        end

    end

    return components;
end

function OAATcclrev( mask::Array{UInt8,3}; connectivity=connectivity8() )

	return OAATccl( mask, connectivity=connectivity, ref=UInt8(0) ) 
end


""" Utility functions """

function CCLlimits( ccl ) 
	miny = typemax(eltype(ccl[1]))
	minx = typemax(eltype(ccl[1]))
	maxy = 0
	maxx = 0
	@simd for coords in ccl
		miny = ( coords[1] < miny ) ? coords[1] : miny
		minx = ( coords[2] < minx ) ? coords[2] : minx
		maxy = ( coords[1] > maxy ) ? coords[1] : maxy
		maxx = ( coords[2] > maxx ) ? coords[2] : maxx
	end
	return miny, minx, maxy, maxx
end

function infoAboutPositiveCCLS( input::Array{<:Real,N}; areaTH=0 ) where {N}

	mean  = 0.0
	@inbounds @simd for idx in 1:length(input)
		mean += input[idx]
	end
	mean  = mean/length(input)
	signs = Array{UInt8,N}(undef,size(input))
	@inbounds @simd for idx in 1:length(input) 
		signs[idx] = UInt8( sign(input[idx]-mean) > 0 )
	end
	
	CCLS    = OAATccl( signs )
    cont    = 0
	maxArea = 0
    for ccl in CCLS
        cont   += length(ccl) > areaTH
		maxArea = ( length(ccl) ) > maxArea ? length(ccl) : maxArea
    end
	return cont, maxArea
end

function countNegativeCCLS( input::Array{<:Real,N}; areaTH=0 ) where {N}

	mean  = 0.0
	@inbounds @simd for idx in 1:length(input)
		mean += input[idx]
	end
	mean = mean/length(input)
	signs = Array{UInt8,N}(undef,size(input))
	@inbounds @simd for idx in 1:length(input) 
		signs[idx] = UInt8( sign(input[idx]-mean) < 0 )
	end
	
	CCLS = OAATccl( signs )
    cont = 0
    for ccl in CCLS
        cont += length(ccl) > areaTH
    end
	return cont
end




# Extra
function OAATcclID( mask::Array{UInt16,N}, ID; connectivity=default_connect(N) ) where {N}
	isID = Array{UInt8,N}(undef,size(mask))
	@inbounds @simd for idx in 1:length(mask)
		isID[idx] == ( mask[idx] == ID )
	end
	return OAATccl( mask, connectivity=connectivity )
end

function OAATccl( mask::Array{UInt16,3};
                  connectivity=( (-1,0,0), (1,0,0), 
                                 (0,-1,0), (0,1,0),
                                 (0,0,-1), (0,0,1) ) )

    h, w, d = size(mask);
    mask_   = copy(mask);
	num     = length( unique( mask[:] ) ); 

    components = [ Array{Tuple{Int64,Int64,Int64},1}(undef,0) for x in 1:num-1 ]
	
	for zet in 1:d, col in 1:w, row in 1:h

        if mask_[ row, col, zet ] > 0

            push!( components[ mask_[row,col,zet] ], (row,col,zet) )
			mask_[row,col,zet] = 0
        end
	end

    return components;
end
