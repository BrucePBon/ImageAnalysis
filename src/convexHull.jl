


function convexHull( points::Array{Tuple{T,T,T},1} ) where {T<:Real}

    @assert length(points) >= 4
    # length(points) for now, to have a static size

    # 3 vertices of each face.
    # The 3 edges of the face are given by faces[idx][1:2], faces[idx][2:3], faces[idx][3]:faces[idx][1];
    faces = Array{Tuple{UInt16,UInt16,UInt16},1}(undef,length(points));
    # Info about neighbouring faces of each face
    neigh = Array{Tuple{UInt16,UInt16,UInt16},1}(undef,length(points));
    # 3 3D normals for each edge of a face, 9 numbers in total
    norms = Array{NTuple{9,Float32},1}(undef,length(points));

    # Intializing first 4 faces
    faces[1] = ( 1, 2, 3 ); neigh[1] = ( 2, 3, 4 );
    faces[2] = ( 2, 3, 4 ); neigh[2] = ( 1, 3, 4 );
    faces[3] = ( 3, 4, 1 ); neigh[3] = ( 1, 2, 4 );
    faces[4] = ( 4, 1, 2 ); neigh[4] = ( 1, 2, 3 );

    # Faces will be modified or added. When adding new faces, "last" keeps track of the index of the last face.
    last = 4;

    center = sum(points[1:4])/4;

    for fc = 1:4
        face = faces[fc]
        # edge between e1 and e2
        for (e1, e2) in [ (1,2), (2,3), (3,1) ]
            p1 = points[ face[e1] ];
            p2 = points[ face[e2] ];
            norms[e1] = normalCenter( p1, p2, center )
        end
    end

    return faces, neigh
end

function normalCenter( p1, p2, ctr )
    v1 =  p2 .- p1;
    v2 = ctr .- p1;
    n  = ctr .- ( v1 .* v2 )/sqrt( sum(v1) )
end
