"""
    FLEXIBLE (BUT NOT NECESSARILY MOST EFFICIENT) METHODS FOR IMPLEMENTING A VARIETY OF MEAN INTENSITY OPERATIONS, 
    SUCH AS MEAN DERIVATIVE COMPUTATIONS. 
"""

# Sampling the mean of the direct neighbouring squares/cubes. 
# Since each neighbouring region will be paired with its opposite (ex (-1,0) with (1,0) ), we only need to 
# define half of the offsets here. The other half can be compute by multiplying by -1. 
direct_offs( ::Val{2} ) = ( (-1,0), (0,-1) ), 2; 
direct_offs( ::Val{3} ) = ( (-1,0,0), (0,-1,0), (0,0,-1) ), 3; 



# Sampling the mean of the direct + diagonal neighboring squares/cubes. Because of pairing of opposite neighbors
# (see above), only half of the neighbours need to be sampled. 
directPlus_offs( ::Val{2} ) = ( (-1,-1), (-1,0), (-1,1), (0,-1) ), 4;  # 8
directPlus_offs( ::Val{3} ) = ( (-1,-1,-1), (-1,0,-1), (-1,1,-1),      # 26
                                ( 0,-1,-1), ( 0,0,-1), ( 0,1,-1),
                                ( 1,-1,-1), ( 1,0,-1), ( 1,1,-1),
                                (-1,-1, 0), (-1,0, 0), (-1,1, 0),
                                ( 0,-1, 0) ), 13; 

# Used in line 39 to figure out the size of the array of offsets in allocating functions.
direct_len( ::Val{2} ) = 4;
direct_len( ::Val{3} ) = 6;
directPlus_len( ::Val{2} ) = 8;
directPlus_len( ::Val{3} ) = 26;

# This function align the pairs of neighbouring regions. Each region's coordinates are followed by 
# their opposite region. This makes it easy to look for maxima or minima, which are characterized by 
# symmetric gradients of intensities, ex grad[2*x] == -1 * grad[2*x-1]. 
# This function also applied a spacing, allowing to scale up the distance at which neighboring 
# means are sampled. 
function get_offsets!( spacing::NTuple{N,<:Integer}, offsets; fun=(x)->direct_offs(x) ) where{ N }
    offs, len = fun( Val(N) );
    for n in 1:len
        offsets[   n*2   ] = offs[n] .* spacing
        offsets[ n*2 - 1 ] = offsets[ n*2 ] .* -1
    end
    return offsets
end

# in-place
directNeighbors_offsets!(     spacing::NTuple{N,<:Integer}, offsets ) where{ N } = get_offsets!( spacing, offsets, fun=(x)->direct_offs(x)   ); 
directPlusNeighbors_offsets!( spacing::NTuple{N,<:Integer}, offsets ) where{ N } = get_offsets!( spacing, offsets, fun=(x)->directPlus_offs(x) ); 

# allocating
directNeighbors_offsets(     spacing::NTuple{N,<:Integer} ) where{ N } = directNeighbors_offsets!( spacing, Array{NTuple{N,Int64},1}(undef,direct_len(Val(N))) );
directPlusNeighbors_offsets( spacing::NTuple{N,<:Integer} ) where{ N } = directPlusNeighbors_offsets!( spacing, Array{NTuple{N,Int64},1}(undef,directPlus_len(Val(N))) ); 

# Apply the offsets to the corners defined by tlf and bre. 
function applyOffsets( tlf, bre, offsets )
    return applyOffsets!( tlf, bre, offsets, [ ( tlf, bre ) for off in offsets ] ); 
end

function applyOffsets!( tlf, bre, offsets, output )
    @inbounds for idx in 1:length(offsets)
        off = offsets[idx]; 
        output[idx] = ( tlf .+ off, bre .+ off )
    end
    return output
end

# After applying the offsets and sampling the sum of values in each neighbouring region from the input integral array
function sample_integralArray( grid_coords, intArray; den=1 )
    results = Array{eltype(intArray),1}(undef,length(grid_coords));
    sample_integralArray!( grid_coords, results, intArray, den=den ); 
end

function sample_integralArray!( grid_coords, results, intArray; den=1 )
    # coords contains the two corners that define the square/cube to sample the sums from. In 2D coords == ( topleft, botright ), and in 3D coords == ( topleftfront, botrightend );
    for gidx in 1:length(grid_coords)
        # Using integralArea_checked to avoid the "headache" of padding the integral array
        results[gidx] = ImageAnalysis.integralArea_checked( intArray, grid_coords[gidx][1], grid_coords[gidx][2] )/den; 
    end
    return results
end

# From the results of sample_integralArrays, computing a score for maxima and minima
countMaxima( sums, ref; f=1 ) = sum( [ (sign(sums[2*n]*f - ref) == -1) && (sign(sums[2*n-1]*f - ref) == -1) for n in 1:div(length(sums),2) ] )
countMinima( sums, ref; f=1 ) = sum( [ (sign(sums[2*n]*f - ref) ==  1) && (sign(sums[2*n-1]*f - ref) ==  1) for n in 1:div(length(sums),2) ] )

function count2ndMaxima( sums1st, sums2nd, ref; f=1 )
    out = 0
    for n in 1:div(length(sums1st),2)
        out += sign( sums1st[2*n]*f - ref ) == -1  && sign( sums2nd[2*n]*f - sums1st[2*n] ) == -1  &&  sign( sums1st[2*n - 1]*f - ref ) == -1  && sign( sums2nd[2*n - 1]*f - sums1st[2*n - 1] ) == -1
    end
    return out
end

function count2ndMinima( sums1st, sums2nd, ref; f=1 )
    out = 0
    for n in 1:div(length(sums1st),2)
        out += ( sign( sums1st[2*n]*f - ref ) == 1  && sign( sums2nd[2*n]*f - sums1st[2*n] ) == 1 )  && ( sign( sums1st[2*n - 1]*f - ref ) == 1 && sign( sums2nd[2*n - 1]*f - sums1st[2*n - 1] ) == 1 )
    end
    return out
end



# Default spacing is the minimum non-overlapping distance, aka "rad2 + rad1 + 1";  
function proto_mean_maxima_pad( input::Array{<:Real,2}, rad1; rad2=rad1, spacing=nothing, threshold=-1, fmax=1 )

    maxima = zeros( Bool, size(input) ); 
    intA   = ImageAnalysis.integralArray( input ); 
    # spacing = offsets of the direct neighbor regions from the central region, allowing to bring them closer or further away from the center.      
    spacing = ( spacing == nothing ) ? rad2 .+ 1 .+ rad1 : spacing; 
    # Array of tuples, with the offsets for each neighbouring region
    # Making a copy, out_offsets, that we can operate in place on
    grid_offsets = directNeighbors_offsets( spacing );
    neigh_coords = applyOffsets( (0,0), (0,0), grid_offsets ); 
    # computing out-of-place once, so we can compute in-place later on. 
    grid_averages = sample_integralArray( neigh_coords, intA ); 
    # we need to divide the sums by the number of pixels to get the average
    npix1 = prod( 2 .* rad1 .+ 1 ); 
    npix2 = prod( 2 .* rad2 .+ 1 ); 
    # for each coordinate of the input, sample its neighbours decide whether the coordinate is the center of an extremum
    for x in 1:size(input,2), y in 1:size(input,1)
        tlf1 = ( y,x ) .- rad1 .- 1; 
        bre1 = ( y,x ) .+ rad1; 
        avg1 = ImageAnalysis.integralArea_checked( intA, tlf1, bre1 )/npix1;   
        ( avg1 <= threshold ) && ( continue; )
        tlf2 = ( y,x ) .- rad2 .- 1; 
        bre2 = ( y,x ) .+ rad2; 
        # applying the neighbouring offsets around the current coordinate
        applyOffsets!( tlf2, bre2, grid_offsets, neigh_coords ); 
        # sampling the sums on the grid from the integral array
        sample_integralArray!( neigh_coords, grid_averages, intA, den=npix2 ); 
        maxscore = countMaxima( grid_averages, avg1 ); 
        minscore = countMinima( grid_averages, avg1 );
        maxima[ y, x ] = maxscore > fmax;
    end
    return maxima
end

# Default spacing is the minimum non-overlapping distance, aka "rad2 + rad1 + 1";  
function proto_mean_maxima_pad( input::Array{<:Real,3}, rad1; rad2=rad1, spacing=nothing, threshold=-1, fmax=1, f=1 )
    maxima = zeros( Bool, size(input) ); 
    intA   = ImageAnalysis.integralArray( input ); 
    spacing = ( spacing == nothing ) ? rad2 .+ 1 .+ rad1 : spacing; 
    grid_offsets = directPlusNeighbors_offsets( spacing ); 
    neigh_coords = applyOffsets( (0,0,0), (0,0,0), grid_offsets ); 
    grid_averages = sample_integralArray( neigh_coords, intA ); 
    npix1 = prod( 2 .* rad1 .+ 1 ); 
    npix2 = prod( 2 .* rad2 .+ 1 ); 
    for z in 1:size(input,3), x in 1:size(input,2), y in 1:size(input,1)
        tlf1 = ( y,x,z ) .- rad1 .- 1; 
        bre1 = ( y,x,z ) .+ rad1; 
        avg1 = ImageAnalysis.integralArea_checked( intA, tlf1, bre1 )/npix1;   
        ( avg1 <= threshold ) && ( continue; )
        tlf2 = ( y,x,z ) .- rad2 .- 1; 
        bre2 = ( y,x,z ) .+ rad2; 
        applyOffsets!( tlf2, bre2, grid_offsets, neigh_coords ); 
        sample_integralArray!( neigh_coords, grid_averages, intA, den=npix2 ); 
        maxscore = countMaxima( grid_averages, avg1, f=f ); 
        maxima[ y, x, z ] = maxscore > fmax;
    end
    return maxima
end

# Default spacing is the minimum non-overlapping distance, aka "rad2 + rad1 + 1";  
function proto_mean_2ndmaxima_pad( input::Array{<:Real,3}, rad1; rad2=rad1, spacing1=nothing, spacing2=nothing, threshold=-1, fmax=1, f=1 )

    maxima = zeros( Bool, size(input) ); 
    intA   = ImageAnalysis.integralArray( input ); 

    spacing1st       = ( spacing1 == nothing ) ? rad2 .+ 1 .+ rad1 : spacing1; 
    grid_offsets1st  = directPlusNeighbors_offsets( spacing1st ); 
    neigh_coords1st  = applyOffsets( (0,0,0), (0,0,0), grid_offsets1st ); 
    grid_averages1st = sample_integralArray( neigh_coords1st, intA ); 

    spacing2nd       = ( spacing2 == nothing ) ? rad2 .+ 1 .+ rad1 : spacing2; 
    grid_offsets2nd  = directPlusNeighbors_offsets( spacing2nd ); 
    neigh_coords2nd  = applyOffsets( (0,0,0), (0,0,0), grid_offsets2nd ); 
    grid_averages2nd = sample_integralArray( neigh_coords2nd, intA ); 

    npix1 = prod( 2 .* rad1 .+ 1 ); 
    npix2 = prod( 2 .* rad2 .+ 1 ); 

    for z in 1:size(input,3), x in 1:size(input,2), y in 1:size(input,1)
        tlf1 = ( y,x,z ) .- rad1 .- 1; 
        bre1 = ( y,x,z ) .+ rad1; 
        avg1 = ImageAnalysis.integralArea_checked( intA, tlf1, bre1 )/npix1;   
        ( avg1 <= threshold ) && ( continue; )
        tlf2 = ( y,x,z ) .- rad2 .- 1; 
        bre2 = ( y,x,z ) .+ rad2; 
        applyOffsets!( tlf2, bre2, grid_offsets1st, neigh_coords1st );
        applyOffsets!( tlf2, bre2, grid_offsets2nd, neigh_coords2nd ); 
        sample_integralArray!( neigh_coords1st, grid_averages1st, intA, den=npix2 ); 
        sample_integralArray!( neigh_coords2nd, grid_averages2nd, intA, den=npix2 ); 
        maxscore = count2ndMaxima( grid_averages1st, grid_averages2nd, avg1, f=f ); 
        maxima[ y, x, z ] = maxscore > fmax;
    end
    return maxima
end

# Default spacing is the minimum non-overlapping distance, aka "rad2 + rad1 + 1";  
function proto_mean_minima_pad( input::Array{<:Real,3}, rad1; rad2=rad1, spacing=nothing, threshold=-1, fmin=1, f=1 )
    minima = zeros( Bool, size(input) ); 
    intA   = ImageAnalysis.integralArray( input ); 
    spacing = ( spacing == nothing ) ? rad2 .+ 1 .+ rad1 : spacing; 
    grid_offsets = directPlusNeighbors_offsets( spacing ); 
    neigh_coords = applyOffsets( (0,0,0), (0,0,0), grid_offsets ); 
    grid_averages = sample_integralArray( neigh_coords, intA ); 
    npix1 = prod( 2 .* rad1 .+ 1 ); 
    npix2 = prod( 2 .* rad2 .+ 1 ); 
    for z in 1:size(input,3), x in 1:size(input,2), y in 1:size(input,1)
        tlf1 = ( y,x,z ) .- rad1 .- 1; 
        bre1 = ( y,x,z ) .+ rad1; 
        avg1 = ImageAnalysis.integralArea_checked( intA, tlf1, bre1 )/npix1;   
        ( avg1 <= threshold ) && ( continue; )
        tlf2 = ( y,x,z ) .- rad2 .- 1; 
        bre2 = ( y,x,z ) .+ rad2; 
        applyOffsets!( tlf2, bre2, grid_offsets, neigh_coords ); 
        sample_integralArray!( neigh_coords, grid_averages, intA, den=npix2 ); 
        maxscore = countMinima( grid_averages, avg1, f=f ); 
        minima[ y, x, z ] = maxscore > fmin;
    end
    return minima
end

# Default spacing is the minimum non-overlapping distance, aka "rad2 + rad1 + 1";  
function proto_mean_minima_pad_masked( input::Array{<:Real,3}, mask, rad1; rad2=rad1, nrad=(2,2,2), spacing=nothing, fmin=1, f=1, volf=0.5 )
    minima = zeros( Bool, size(input) ); 
    intA   = ImageAnalysis.integralArray( input );
    intN   = ImageAnalysis.integralArray( mask ); 
    spacing = ( spacing == nothing ) ? rad2 .+ 1 .+ rad1 : spacing; 
    grid_offsets = directPlusNeighbors_offsets( spacing ); 
    neigh_coords = applyOffsets( (0,0,0), (0,0,0), grid_offsets ); 
    grid_averages = sample_integralArray( neigh_coords, intA ); 
    npix1 = prod( 2 .* rad1 .+ 1 ); 
    npix2 = prod( 2 .* rad2 .+ 1 ); 
    for z in 1:size(input,3), x in 1:size(input,2), y in 1:size(input,1)
        if mask[y,x,z]
            tlf1 = ( y,x,z ) .- rad1 .- 1; 
            bre1 = ( y,x,z ) .+ rad1; 
            avg1 = ImageAnalysis.integralArea_checked( intA, tlf1, bre1 )/npix1;   
            tlfn = ( y,x,z ) .- nrad .- 1; 
            bren = ( y,x,z ) .+ nrad; 
            N1   = ImageAnalysis.integralArea_checked( intN, tlfn, bren );   
            tlf2 = ( y,x,z ) .- rad2 .- 1; 
            bre2 = ( y,x,z ) .+ rad2; 
            applyOffsets!( tlf2, bre2, grid_offsets, neigh_coords ); 
            sample_integralArray!( neigh_coords, grid_averages, intA, den=npix2 ); 
            minscore = countMinima( grid_averages, avg1, f=f ); 
            minima[ y, x, z ] = minscore > fmin && N1 < volf;
        end
    end
    return minima
end

# Default spacing is the minimum non-overlapping distance, aka "rad2 + rad1 + 1";  
function proto_mean_2ndminima_pad_masked( input::Array{<:Real,3}, mask, rad1; rad2=rad1, spacing1=nothing, spacing2=nothing, fmin=1, f=1 )
    minima = zeros( Bool, size(input) ); 
    intA   = ImageAnalysis.integralArray( input ); 
    intN   = ImageAnalysis.integralArray( mask  ); 

    spacing1st       = ( spacing1 == nothing ) ? rad2 .+ 1 .+ rad1 : spacing1; 
    grid_offsets1st  = directPlusNeighbors_offsets( spacing1st ); 
    neigh_coords1st  = applyOffsets( (0,0,0), (0,0,0), grid_offsets1st ); 
    grid_averages1st = sample_integralArray( neigh_coords1st, intA ); 
    grid_N1st        = sample_integralArray( neigh_coords1st, intN ); 


    spacing2nd       = ( spacing2 == nothing ) ? rad2 .+ 1 .+ rad1 : spacing2; 
    grid_offsets2nd  = directPlusNeighbors_offsets( spacing2nd ); 
    neigh_coords2nd  = applyOffsets( (0,0,0), (0,0,0), grid_offsets2nd ); 
    grid_averages2nd = sample_integralArray( neigh_coords2nd, intA ); 
    grid_N2nd        = sample_integralArray( neigh_coords2nd, intN ); 

    npix1 = prod( 2 .* rad1 .+ 1 ); 
    npix2 = prod( 2 .* rad2 .+ 1 ); 

    for z in 1:size(input,3), x in 1:size(input,2), y in 1:size(input,1)
        if mask[y,x,z]
            tlf1 = ( y,x,z ) .- rad1 .- 1; 
            bre1 = ( y,x,z ) .+ rad1; 
            avg1 = ImageAnalysis.integralArea_checked( intA, tlf1, bre1 )/npix1;#/ImageAnalysis.integralArea_checked( intN, tlf1, bre1 );   

            tlf2 = ( y,x,z ) .- rad2 .- 1; 
            bre2 = ( y,x,z ) .+ rad2; 
            applyOffsets!( tlf2, bre2, grid_offsets1st, neigh_coords1st );
            applyOffsets!( tlf2, bre2, grid_offsets2nd, neigh_coords2nd ); 
            sample_integralArray!( neigh_coords1st, grid_averages1st, intA, den=npix2 ); 
            sample_integralArray!( neigh_coords1st, grid_N1st, intN ); 
            sample_integralArray!( neigh_coords2nd, grid_averages2nd, intA, den=npix2 );
            sample_integralArray!( neigh_coords2nd, grid_N2nd, intN );

            #grid_averages1st ./= grid_N1st
            #grid_averages2nd ./= grid_N2nd
 
            minscore = count2ndMinima( grid_averages1st, grid_averages2nd, avg1, f=f ); 
            minima[ y, x, z ] = minscore > fmin;
        end
    end
    return minima
end

# Default spacing is the minimum non-overlapping distance, aka "rad2 + rad1 + 1";  
function proto_mean_2ndminima_pad_masked_iso( input::Array{<:Real,3}, mask, rad1; rad2=rad1.*2, rad3=rad1.*3, fmin=1, f1=1, f2=1 )

    minima = zeros( Bool, size(input) ); 
    intA   = ImageAnalysis.integralArray( input ); 
    intN   = ImageAnalysis.integralArray( mask  ); 


    for z in 1:size(input,3), x in 1:size(input,2), y in 1:size(input,1)
        if mask[y,x,z]
            tlf1 = ( y,x,z ) .- rad1 .- 1; 
            bre1 = ( y,x,z ) .+ rad1; 
            sum1 = ImageAnalysis.integralArea_checked( intA, tlf1, bre1 );
            N1   = ImageAnalysis.integralArea_checked( intN, tlf1, bre1 ); 
            avg1 = sum1 / N1;
            
            tlf2 = ( y,x,z ) .- rad2 .- 1; 
            bre2 = ( y,x,z ) .+ rad2; 
            sum2 = ImageAnalysis.integralArea_checked( intA, tlf2, bre2 ) - sum1;
            N2   = ImageAnalysis.integralArea_checked( intN, tlf2, bre2 ) - N1;   
            avg2 = sum2 / N2; 

            tlf3 = ( y,x,z ) .- rad3 .- 1; 
            bre3 = ( y,x,z ) .+ rad3; 
            sum3 = ImageAnalysis.integralArea_checked( intA, tlf3, bre3 ) - sum2 - sum1; 
            N3   = ImageAnalysis.integralArea_checked( intN, tlf3, bre3 ) - N2 - N1;   
            avg3 = sum3 / N3; 
 
            minima[ y, x, z ] = ( avg3 > avg2*f2 ) + ( avg2*f2 > avg1*f1 ) > fmin;
        end
    end
    return minima
end

# Default spacing is the minimum non-overlapping distance, aka "rad2 + rad1 + 1";  
function proto_mean_maxima_pad_iso( input::Array{<:Real,3}, rad1; rad2=rad1.*2, f=1, threshold=0 )

    maxima = zeros( Bool, size(input) ); 
    intA   = ImageAnalysis.integralArray( input ); 

    npix1  = prod( 2 .* rad1 .+ 1 ); 
    npix2  = prod( 2 .* rad2 .+ 1 ); 


    for z in 1:size(input,3), x in 1:size(input,2), y in 1:size(input,1)

        tlf1 = ( y,x,z ) .- rad1 .- 1; 
        bre1 = ( y,x,z ) .+ rad1; 
        sum1 = ImageAnalysis.integralArea_checked( intA, tlf1, bre1 );
        avg1 = sum1 / npix1;

        if avg1 < threshold
            continue
        end
        
        tlf2 = ( y,x,z ) .- rad2 .- 1; 
        bre2 = ( y,x,z ) .+ rad2; 
        sum2 = ImageAnalysis.integralArea_checked( intA, tlf2, bre2 );
        avg2 = sum2 / npix2; 

        maxima[ y, x, z ] = avg1/avg2 > f;
    end
    return maxima
end

# Default spacing is the minimum non-overlapping distance, aka "rad2 + rad1 + 1";  
function proto_mean_maxima_pad_masked_iso( input::Array{<:Real,3}, mask, rad1; rad2=rad1.*2, f=1, threshold=0 )

    maxima = zeros( Bool, size(input) ); 
    intA   = ImageAnalysis.integralArray( input ); 
    intN   = ImageAnalysis.integralArray(  mask ); 
    for z in 1:size(input,3), x in 1:size(input,2), y in 1:size(input,1)
        if mask[y,x,z]
            tlf1 = ( y,x,z ) .- rad1 .- 1; 
            bre1 = ( y,x,z ) .+ rad1; 
            sum1 = ImageAnalysis.integralArea_checked( intA, tlf1, bre1 );
            N1   = ImageAnalysis.integralArea_checked( intN, tlf1, bre1 );
            avg1 = sum1 / N1;

            if avg1 < threshold
                continue
            end
            
            tlf2 = ( y,x,z ) .- rad2 .- 1; 
            bre2 = ( y,x,z ) .+ rad2; 
            sum2 = ImageAnalysis.integralArea_checked( intA, tlf2, bre2 );
            N2   = ImageAnalysis.integralArea_checked( intN, tlf2, bre2 );
            avg2 = sum2 / N2; 

            maxima[ y, x, z ] = avg1/avg2 > f;
        end
    end
    return maxima
end

# Default spacing is the minimum non-overlapping distance, aka "rad2 + rad1 + 1";  

proto_mean_opening( mask, rad; th=0 ) = proto_mean_opening!( mask, ImageAnalysis.integralArray( mask ), rad, th=th ); 

function proto_mean_opening!( mask::Array{<:Integer,3}, intN, rad; th=0 )

    opened = zeros( Bool, size(mask) ); 
    for z in 1:size(mask,3), x in 1:size(mask,2), y in 1:size(mask,1)
       #( !mask[y,x,z] ) && ( continue; )
        N1 = ImageAnalysis.integralArea_checked( intN, ( y,x,z ) .- rad .- 1, ( y,x,z ) .+ rad );
        opened[ y, x, z ] = N1 >= th ;
    end
    return opened
end

proto_any_opening( mask, rad, th=0 ) = proto_any_opening!( ImageAnalysis.integralArray( mask ), rad, th=th ); 

function proto_any_opening!( mask, intN, rad; th=0 )
    for z in 1:size(intN,3)-1, x in 1:size(intN,2)-1, y in 1:size(intN,1)-1
        ( !mask[y,x,z] ) && continue; 
        N1 = ImageAnalysis.integralArea_checked( intN, ( y,x,z ) .- rad .- 1, ( y,x,z ) .+ rad );
        ( N1 > th ) && return true
    end
    return false
end

# TODO: needs rethinking

function proto_mean_opening_iterative( mask::Array{<:Integer,3}, rad; max_th=100, min_th=0, step=5, stop=20)

    intN = ImageAnalysis.integralArray( mask );
    t0 = min_th
    t1 = max_th
    while ( t1 - t0 ) > 1
        tmp = div( t0 + t1, 2 ); 
        if proto_any_opening!( mask, intN, rad; th=tmp )
            t0 = tmp
        else
            t1 = tmp
        end
    end
    #println( " ths: ", proto_any_opening!( mask, intN, rad; th=t0 ), "  ", proto_any_opening!( mask, intN, rad; th=t1 ) )

    println( t0:t1 )

    peaks = proto_mean_opning( mask, rad, t0 ); 

    finalLabels = ImageComponentAnalysis.label_compoments( peaks );
    lastIdx     = maximum( finalLabels ); 
    intNLabels  = zeros( Float32, size( mask ) .+ 1 ); 

    # keeping track of voxels labels and number of labels
    intLabels   = zeros( Float32, size( mask ) .+ 1 ); 
    NLabels     = zeros( Float32, size( mask ) ); 

    for th in t0:t0
        lastIdx = proto_add_count_labels!( mask, intN, intLabels, intNLabels, NLabels, finalLabels, th, rad, lastIdx )
    end

    return finalLabels
end

# TODO: needs rethinking

function proto_add_count_labels!( mask, intMask, intlabels, intNlabels, Nlabels, finallabels, th, rad, lastIdx )

    nnew  = 0; 
    njoin = 0;

    for z in size(intMask,3)-1:-1:1, x in size(intMask,2)-1:-1:1, y in size(intMask,1)-1:-1:1

        ( !mask[y,x,z] ) && continue; 

        Nmask = ImageAnalysis.integralArea_checked( intMask, ( y,x,z ) .- rad .- 1, ( y,x,z ) .+ rad );
        Nlbls = ImageAnalysis.integralArea_checked( intNlabels, ( y,x,z ) .- rad .- 1, ( y,x,z ) .+ rad );

        # if number of neighbours > th
        if ( Nmask > th )
            # if the neighbour is already labelled, take the label of the neighbour
            if ( Nlbls == 1 )
                finallabels[y,x,z] = ImageAnalysis.integralArea_checked( intlabels, ( y,x,z ) .- rad .- 1, ( y,x,z ) .+ rad )
                mask[y,x,z] = false
                njoin += 1
            else
            # it has no neighbours, start a new label
                # intLabels[ 2:y+1, 2:x+1, 2:z+1 ] .+= lastIdx # slow part
                #intNlabels[ 2:y+1, 2:x+1, 2:z+1 ] .+= 1       # slow part
                mask[y,x,z] = false
                finallabels[y,x,z] = lastIdx
                Nlabels[y,x,z] = 1
                lastIdx += 1
                nnew += 1
            end
        end
    end

    println( "nnew, njoin: ", ( nnew, njoin ) )

    integralArray!( Nlabels, intNlabels )
    integralArray!( finallabels, intlabels )

    return lastIdx
end
