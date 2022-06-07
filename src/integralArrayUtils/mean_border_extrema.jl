# What fraction of asterisks (*) have a lower mean that X

# * * * * *      omul = 2
# * x x x *
# * X X X * 
# * * * * *

# osize, size of the "asterisks", which are square/cubic subdivisions
# isize, inner size, is the size of the X-covered region. To keep things
#        easy, isize is a multiple of osize. 

function mean_border_extrema( vol, omul, osize, f )

    # inner size
    isize = osize .* omul; 
    
    # osize | vol | osize
    padsize = size( vol ) .+ 2 .* osize;
    intA    = zeros( Float32, padsize .+ 1 );
    integralArray_pad!( vol, intA, osize  ); 
    
    # convenient quantities
    lows  = (  1, 1, 1 ) .+ osize .+ isize; 
    highs =  size( vol ) .+ osize;

    # outer sums (an inner). TODO: optimize
    odims = prod( omul .+ 2 ) .- prod( omul )
    osums = zeros( Float32, odims ); 
    score = zeros( Float32, size( vol ) ); 

    for z in lows[3]:highs[3], x in lows[2]:highs[2], y in lows[1]:highs[1]

        itl  = ( y, x, z ) .- isize .- 1; 
        ibr  = ( y, x, z ); 
        isum = ImageAnalysis.integralArea( intA, itl, ibr ) / prod( omul );  
        
        cont = 1; 
        # 1st
        for xstep in -1:omul[2], ystep in -1:omul[1]
            otl  = itl .+ ( ystep, xstep, -1 ) .* osize
            obr  = otl .+ osize .- 1; 
            osum = ImageAnalysis.integralArea( intA, otl, obr ); 
            osums[cont] = osum;   
            cont += 1
        end
        # 2nd
        for zstep in 0:omul[3]-1, xstep in -1:omul[2]
            yrange = ( xstep == -1 || xstep == omul[2] ) ?  UnitRange(-1,omul[1]) : [-1,omul[1]];
            for ystep in yrange
                otl  = itl .+ ( ystep, xstep, zstep ) .* osize
                obr  = otl .+ osize .- 1
                osum = ImageAnalysis.integralArea( intA, otl, obr ); 
                osums[cont] = osum;   
                cont += 1
            end
        end
        # 3rd
        for xstep in -1:omul[2], ystep in -1:omul[1]
            otl  = itl .+ ( ystep, xstep, omul[3] ) .* osize
            obr  = otl .+ osize .- 1
            osum = ImageAnalysis.integralArea( intA, otl, obr ); 
            osums[cont] = osum;   
            cont += 1
        end 

        hits = 0
        for idx in 1:length(osums)
            hits += osums[idx] < isum
        end

        score[ y-osize[1]-isize[1], x-osize[2]-isize[2], z-osize[3]-isize[3] ] = hits
    end

    return score
end


function mean_sums( vol, ssize )

    # inner size
    isize = osize .* omul; 
    
    # osize | vol | osize
    padsize = size( vol ) .+ 2 .* osize;
    intA    = zeros( Float32, padsize .+ 1 );
    integralArray_pad!( vol, intA, osize  ); 
    
    # convenient quantities
    lows  = (  1, 1, 1 ) .+ osize .+ isize; 
    highs =  size( vol ) .+ osize;

    # outer sums (an inner). TODO: optimize
    odims = prod( omul .+ 2 ) .- prod( omul )
    osums = zeros( Float32, odims ); 
    score = zeros( Float32, size( vol ) ); 

    for z in lows[3]:highs[3], x in lows[2]:highs[2], y in lows[1]:highs[1]

        itl  = ( y, x, z ) .- isize .- 1; 
        ibr  = ( y, x, z ); 
        isum = ImageAnalysis.integralArea( intA, itl, ibr ) / prod( omul );  
        
        cont = 1; 
        # 1st
        for xstep in -1:omul[2], ystep in -1:omul[1]
            otl  = itl .+ ( ystep, xstep, -1 ) .* osize
            obr  = otl .+ osize .- 1; 
            osums[cont] = ImageAnalysis.integralArea( intA, otl, obr ); 
            cont += 1
        end
        # 2nd
        for zstep in 0:omul[3]-1, xstep in -1:omul[2]
            yrange = ( xstep == -1 || xstep == omul[2] ) ?  UnitRange(-1,omul[1]) : [-1,omul[1]];
            for ystep in yrange
                otl  = itl .+ ( ystep, xstep, zstep ) .* osize
                obr  = otl .+ osize .- 1
                osums[cont] = ImageAnalysis.integralArea( intA, otl, obr ); 
                cont += 1
            end
        end
        # 3rd
        for xstep in -1:omul[2], ystep in -1:omul[1]
            otl  = itl .+ ( ystep, xstep, omul[3] ) .* osize
            obr  = otl .+ osize .- 1
            osums[cont] = ImageAnalysis.integralArea( intA, otl, obr );    
            cont += 1
        end 

        hits = 0
        for idx in 1:length(osums)
            hits += osums[idx] < isum
        end

        score[ y-osize[1]-isize[1], x-osize[2]-isize[2], z-osize[3]-isize[3] ] = hits
    end

    return score
end

