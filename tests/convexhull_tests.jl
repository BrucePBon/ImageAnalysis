using Revise
using ImageAnalysis
using Test

""" testing function to determine if two points are equal """
function test_areEqual( ; N=1, dims=2 ) 
    points  = [ Tuple(rand(dims)) for n in 1:N ]; 
    results = [ ImageAnalysis.areEqual( p, p ) for p in points ]
    return all( results )
end


createColinearPoint( x, y ) = x .+ ( x .- y ) .* rand()

""" testing function to determine if three points are colinear """
function test_areColinear( ; N=1, dims=2 )
    x = [ Tuple(rand(dims)) for n in 1:N ]
    y = [ Tuple(rand(dims)) for n in 1:N ]
    z = createColinearPoint.( x, y ) 
    results = [ ImageAnalysis.areColinear( x[n], y[n], z[n] ) for n in 1:N ]
    return all( results )
end


createCoplanarPoint( x, y, z ) = x .+ ( x .- y ) .* rand() .+ ( x .- z ) .* rand(); 

""" testing function to determine if four points are coplanar """
function test_areCoplanar( ; N=1, dims=2 )
    x = [ Tuple(rand(dims)) for n in 1:N ]
    y = [ Tuple(rand(dims)) for n in 1:N ]
    z = [ Tuple(rand(dims)) for n in 1:N ]
    p = createCoplanarPoint.( x, y, z ) 
    results = [ ImageAnalysis.areCoplanar( x[n], y[n], z[n], p[n] ) for n in 1:N ]
    return all( results )
end


triangle_points( rad ) = [ (rad*sind(90+120),rad*cosd(90+120)), (rad*sind(90),rad*cosd(90)), (rad*sind(90-120),rad*cosd(90-120)) ]; 

""" The 2DConvexHull must include the 3 outer points of a triangle. """
function test_2DConvexHull_empty_triangle(; N=10 )
    points       = triangle_points( 1 ); 
    edges, norms = ImageAnalysis.convexHull( points ); 
    edgepoints   = points[ unique( Base.vcat( [ x[1] for x in edges ], [ x[2] for x in edges ] ) ) ]
    return all( [ p in edgepoints for p in points ] )
end

function convexWeights( N )
    weights = rand( N ); 
    return weights ./ sum( weights )
end

function filled_triangle_points( rad, numpoints ) 
    trig         = Base.vcat( triangle_points( rad ), [(0.0,0.0),] )
    weights      = [ convexWeights( 4 ) for n in 1:numpoints ]
    innerpoints  = [ w[1].*trig[1] .+ w[2].*trig[2] .+ w[3].*trig[3] .+ w[4].*trig[4] for w in weights  ]
    return Base.vcat( innerpoints, trig )
end

""" The 2DConvexHull must include the 3 outer points of a triangle. """
function test_2DConvexHull_filled_triangle(; N=10 )
    trig         = triangle_points( 1 ); 
    points       = filled_triangle_points( 1, 10 );
    edges, norms = ImageAnalysis.convexHull( points ); 
    edgepoints   = points[ unique( Base.vcat( [ x[1] for x in edges ], [ x[2] for x in edges ] ) ) ]
    return all( [ p in edgepoints for p in trig ] );
end

""" .............................................TESTING.................................................. """

println( "TESTING" )
@testset "test" begin

    """ TESTING OF THE GEOMETRIC OPERATIONS INVOLVED IN THE CONVEX HULL ALGORITHM """
    # areSame( x, y ) = ( x == y )
    @test test_areEqual() N=1000 dims=2
    @test test_areEqual() N=1000 dims=3
    # areColinear( x, y, z ) = abs( normdot( x .- z, y .- z )  ) approx 1
    @test test_areColinear() N=100 dims=2
    @test test_areColinear() N=100 dims=3
    # areCoplanar( x, y, z, p ) = dot( x .- p, cross( y .- p, z -. p ) ) == 0
    @test test_areCoplanar() N=100 dims=3

    """ TESTING SIMPLE CASES OF 2D CONVEX HULL """
    # The convexHull of an empty triangle should contain the traingle's 3 corners
    @test test_2DConvexHull_empty_triangle()
    # The ConvexHull of a filled triangle should contain the traingle's 3 corners
    @test test_2DConvexHull_filled_triangle()
end