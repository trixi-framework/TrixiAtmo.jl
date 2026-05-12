using PairedLinkedLists: DoublyLinkedList, getnode
using PairedLinkedLists
using LinearAlgebra

import Base: -
import Base: /
import Base: +
import Base: *
import LinearAlgebra: norm
import LinearAlgebra: dot

mutable struct Point
    x::Float64
    y::Float64
    z::Float64
end

abstract type ElementType end
struct Tri <: ElementType end

abstract type GridForm end
struct CartesianGrid <: GridForm end
struct SphericalGrid <: GridForm end
#=
struct ListNode{T}
    value::T
    next::Union{ListNode{T}, Nothing}
end
=#
mutable struct LinkedList{T}
    head::Union{ListNode{T}, Nothing}
end

function Norm(P::Point)
    sqrt(P.x * P.x + P.y * P.y + P.z * P.z)
end

function Div(P::Point, s)
    return Point(P.x / s, P.y / s, P.z / s)
end

function MidPoint(P1::Point, P2::Point)
    P = Point(0.5 * (P1.x + P2.x),
              0.5 * (P1.y + P2.y), 0.5 * (P1.z + P2.z))
    M = Div(P, Norm(P))
end

function Point(P::Array{Float64, 1})
    return Point(P[1],
                 P[2],
                 P[3])
end
-(P::Point) = Point([-P.x, -P.y, -P.z])
-(P1::Point, P2::Point) = Point([P1.x - P2.x, P1.y - P2.y, P1.z - P2.z])
+(P1::Point, P2::Point) = Point([P1.x + P2.x, P1.y + P2.y, P1.z + P2.z])
/(P::Point, s::Float64) = Point([P.x / s, P.y / s, P.z / s])
*(s::Float64, P::Point) = Point([s * P.x, s * P.y, s * P.z])
*(P::Point, s::Float64) = Point([s * P.x, s * P.y, s * P.z])
*(s::Float32, P::Point) = Point([s * P.x, s * P.y, s * P.z])
*(P::Point, s::Float32) = Point([s * P.x, s * P.y, s * P.z])
norm(P::Point) = sqrt(P.x * P.x + P.y * P.y + P.z * P.z)
dot(P1::Point, P2::Point) = P1.x * P2.x + P1.y * P2.y + P1.z * P2.z
function cross(P1::Point, P2::Point)
    Point([P1.y * P2.z - P1.z * P2.y, P1.z * P2.x - P1.x * P2.z, P1.x * P2.y - P1.y * P2.x])
end

function normalize(P)
    P = P / norm(P)
end

mutable struct NodeTri_T
    P::Point
    Edge::Array{Int, 1}
    Number::Int
end

function NodeTri_T(P)
    Edge = zeros(Int, 0)
    NodeTri_T(P, Edge, 0)
end

mutable struct EdgeTri_T{TN}
    Node1::ListNode{TN}
    Node2::ListNode{TN}
    Face::Array{Int, 1}
    Number::Int
end

function EdgeTri_T(Node1, Node2)
    Face = zeros(Int, 0)
    EdgeTri_T(Node1, Node2, Face, 0)
end

mutable struct FaceTri_T{TE}
    Edge1::ListNode{TE}
    Edge2::ListNode{TE}
    Edge3::ListNode{TE}
    OrientE1::Int
    OrientE2::Int
    OrientE3::Int
    Number::Int
end

function FaceTri_T(Edge1, Edge2, Edge3, OrientE1, OrientE2, OrientE3)
    FaceTri_T(Edge1, Edge2, Edge3, OrientE1, OrientE2, OrientE3, 0)
end

mutable struct TriangularGrid_T{TN, TE, TF}
    NodeList::DoublyLinkedList{TN}
    EdgeList::DoublyLinkedList{TE}
    FaceList::DoublyLinkedList{TF}
    NumNodes::Int
    NumEdges::Int
    NumFaces::Int
end

function TriangularGrid_T()
    NodeList = DoublyLinkedList{NodeTri_T}()
    EdgeList = DoublyLinkedList{EdgeTri_T}()
    FaceList = DoublyLinkedList{FaceTri_T}()
    NumNodes = 0
    NumEdges = 0
    NumFaces = 0
    return TriangularGrid_T{NodeTri_T, EdgeTri_T, FaceTri_T}(NodeList,
                                                             EdgeList,
                                                             FaceList,
                                                             NumNodes,
                                                             NumEdges,
                                                             NumFaces)
end

function CreateIcosahedronGrid()
    Rad = 1.0
    IcosahedronGrid = TriangularGrid_T()

    # Nodes
    Node = NodeTri_T(Point(0.0, 0.0, Rad))
    push!(IcosahedronGrid.NodeList, Node)
    phi = atan(0.5)
    lam = 0.0
    for i in 1:5
        Node = NodeTri_T(Point(Rad * cos(lam) * cos(phi), Rad * sin(lam) * cos(phi),
                               Rad * sin(phi)))
        push!(IcosahedronGrid.NodeList, Node)
        lam = lam + 2.0 * pi / 5.0
    end
    phi = -atan(0.5)
    lam = pi / 5.0
    for i in 1:5
        Node = NodeTri_T(Point(Rad * cos(lam) * cos(phi), Rad * sin(lam) * cos(phi),
                               Rad * sin(phi)))
        push!(IcosahedronGrid.NodeList, Node)
        lam = lam + 2.0 * pi / 5.0
    end
    Node = NodeTri_T(Point(0.0, 0.0, -Rad))
    push!(IcosahedronGrid.NodeList, Node)

    NodeTop = getnode(IcosahedronGrid.NodeList, 1)
    NodeLayerTopFirst = getnode(IcosahedronGrid.NodeList, 2)
    NodeLayerBottomFirst = getnode(IcosahedronGrid.NodeList, 7)
    NodeBottom = getnode(IcosahedronGrid.NodeList, 12)

    # Edges
    NodeLayerTop = NodeLayerTopFirst
    for i in 1:5
        Edge = EdgeTri_T(NodeLayerTop, NodeTop)
        push!(IcosahedronGrid.EdgeList, Edge)
        NodeLayerTop = NodeLayerTop.next
    end
    NodeLayerTop = NodeLayerTopFirst
    for i in 1:4
        Edge = EdgeTri_T(NodeLayerTop, NodeLayerTop.next)
        push!(IcosahedronGrid.EdgeList, Edge)
        NodeLayerTop = NodeLayerTop.next
    end
    Edge = EdgeTri_T(NodeLayerTop, NodeLayerTopFirst)
    push!(IcosahedronGrid.EdgeList, Edge)

    NodeLayerTop = NodeLayerTopFirst
    NodeLayerBottom = NodeLayerBottomFirst
    for i in 1:4
        Edge = EdgeTri_T(NodeLayerBottom, NodeLayerTop)
        push!(IcosahedronGrid.EdgeList, Edge)
        Edge = EdgeTri_T(NodeLayerBottom, NodeLayerTop.next)
        push!(IcosahedronGrid.EdgeList, Edge)
        NodeLayerTop = NodeLayerTop.next
        NodeLayerBottom = NodeLayerBottom.next
    end
    Edge = EdgeTri_T(NodeLayerBottom, NodeLayerTop)
    push!(IcosahedronGrid.EdgeList, Edge)
    Edge = EdgeTri_T(NodeLayerBottom, NodeLayerTopFirst)
    push!(IcosahedronGrid.EdgeList, Edge)

    NodeLayerBottom = NodeLayerBottomFirst
    for i in 1:4
        Edge = EdgeTri_T(NodeLayerBottom, NodeLayerBottom.next)
        push!(IcosahedronGrid.EdgeList, Edge)
        NodeLayerBottom = NodeLayerBottom.next
    end
    Edge = EdgeTri_T(NodeLayerBottom, NodeLayerBottomFirst)
    push!(IcosahedronGrid.EdgeList, Edge)

    NodeLayerBottom = NodeLayerBottomFirst
    for i in 1:5
        Edge = EdgeTri_T(NodeBottom, NodeLayerBottom)
        push!(IcosahedronGrid.EdgeList, Edge)
        NodeLayerBottom = NodeLayerBottom.next
    end

    EdgeTopFirst = getnode(IcosahedronGrid.EdgeList, 1)
    EdgeLayerTopFirst = getnode(IcosahedronGrid.EdgeList, 6)
    EdgeMidFirst = getnode(IcosahedronGrid.EdgeList, 11)
    EdgeLayerBottomFirst = getnode(IcosahedronGrid.EdgeList, 21)
    EdgeBottomFirst = getnode(IcosahedronGrid.EdgeList, 26)

    #Faces
    EdgeTop = EdgeTopFirst
    EdgeLayerTop = EdgeLayerTopFirst
    for i in 1:4
        Face = FaceTri_T(EdgeLayerTop, EdgeTop.next, EdgeTop, 1, 1, -1)
        push!(IcosahedronGrid.FaceList, Face)
        EdgeTop = EdgeTop.next
        EdgeLayerTop = EdgeLayerTop.next
    end
    Face = FaceTri_T(EdgeLayerTop, EdgeTopFirst, EdgeTop, 1, 1, -1)
    push!(IcosahedronGrid.FaceList, Face)

    EdgeLayerTop = EdgeLayerTopFirst
    EdgeMid = EdgeMidFirst
    EdgeLayerBottom = EdgeLayerBottomFirst
    for i in 1:4
        Face = FaceTri_T(EdgeMid, EdgeMid.next, EdgeLayerTop, -1, 1, -1)
        push!(IcosahedronGrid.FaceList, Face)
        EdgeMid = EdgeMid.next
        Face = FaceTri_T(EdgeMid.next, EdgeMid, EdgeLayerBottom, 1, -1, 1)
        push!(IcosahedronGrid.FaceList, Face)
        EdgeLayerTop = EdgeLayerTop.next
        EdgeMid = EdgeMid.next
        EdgeLayerBottom = EdgeLayerBottom.next
    end
    Face = FaceTri_T(EdgeMid, EdgeMid.next, EdgeLayerTop, -1, 1, -1)
    push!(IcosahedronGrid.FaceList, Face)
    EdgeMid = EdgeMid.next
    Face = FaceTri_T(EdgeMidFirst, EdgeMid, EdgeLayerBottom, 1, -1, 1)
    push!(IcosahedronGrid.FaceList, Face)

    EdgeLayerBottom = EdgeLayerBottomFirst
    EdgeBottom = EdgeBottomFirst
    for i in 1:4
        Face = FaceTri_T(EdgeLayerBottom, EdgeBottom.next, EdgeBottom, 1, -1, 1)
        push!(IcosahedronGrid.FaceList, Face)
        EdgeLayerBottom = EdgeLayerBottom.next
        EdgeBottom = EdgeBottom.next
    end
    Face = FaceTri_T(EdgeLayerBottom, EdgeBottomFirst, EdgeBottom, 1, -1, 1)
    push!(IcosahedronGrid.FaceList, Face)
    return IcosahedronGrid
end

function RefineEdge!(Edge, NodeList, EdgeList)
    P1 = Edge.data.Node1.data.P
    P2 = Edge.data.Node2.data.P
    NodeM = newnode(NodeList, NodeTri_T(MidPoint(P1, P2)))
    insertafter!(NodeM, Edge.data.Node1)
    EdgeNew = newnode(EdgeList, EdgeTri_T(NodeM, Edge.data.Node2))
    insertafter!(EdgeNew, Edge)
    Edge.data.Node2 = NodeM
end

function RefineEdgeTriangularGrid!(TriangularGrid)
    NodeList = TriangularGrid.NodeList
    EdgeList = TriangularGrid.EdgeList

    Edge = head(TriangularGrid.EdgeList)
    while ~attail(Edge)
        RefineEdge!(Edge, NodeList, EdgeList)
        Edge = Edge.next.next
    end
end

function RefineFace!(Face, EdgeList, FaceList)
    Edge1 = Face.data.Edge1
    Edge1N = Edge1.next
    Edge2 = Face.data.Edge2
    Edge2N = Edge2.next
    Edge3 = Face.data.Edge3
    Edge3N = Edge3.next
    if Face.data.OrientE1 == 1
        Face2Edge2 = Edge1
        Face2OrientE2 = 1
        Face3Edge1 = Edge1N
        Face3OrientE1 = 1
    else
        Face2Edge2 = Edge1N
        Face2OrientE2 = -1
        Face3Edge1 = Edge1
        Face3OrientE1 = -1
    end
    if Face.data.OrientE2 == 1
        Face3Edge2 = Edge2
        Face3OrientE2 = 1
        Face1Edge1 = Edge2N
        Face1OrientE1 = 1
    else
        Face3Edge2 = Edge2N
        Face3OrientE2 = -1
        Face1Edge1 = Edge2
        Face1OrientE1 = -1
    end
    if Face.data.OrientE3 == 1
        Face1Edge2 = Edge3
        Face1OrientE2 = 1
        Face2Edge1 = Edge3N
        Face2OrientE1 = 1
    else
        Face1Edge2 = Edge3N
        Face1OrientE2 = -1
        Face2Edge1 = Edge3
        Face2OrientE1 = -1
    end
    EdgeI1 = newnode(EdgeList, EdgeTri_T(Edge2.data.Node2, Edge3.data.Node2))
    insertafter!(EdgeI1, Edge1N)
    EdgeI2 = newnode(EdgeList, EdgeTri_T(Edge3.data.Node2, Edge1.data.Node2))
    insertafter!(EdgeI2, Edge2N)
    EdgeI3 = newnode(EdgeList, EdgeTri_T(Edge1.data.Node2, Edge2.data.Node2))
    insertafter!(EdgeI3, Edge3N)

    Face1Edge3 = EdgeI1
    Face1OrientE3 = -1
    Face1 = newnode(FaceList,
                    FaceTri_T(Face1Edge1, Face1Edge2, Face1Edge3, Face1OrientE1,
                              Face1OrientE2, Face1OrientE3))
    insertafter!(Face1, Face)

    Face2Edge3 = EdgeI2
    Face2OrientE3 = -1
    Face2 = newnode(FaceList,
                    FaceTri_T(Face2Edge1, Face2Edge2, Face2Edge3, Face2OrientE1,
                              Face2OrientE2, Face2OrientE3))
    insertafter!(Face2, Face)

    Face3Edge3 = EdgeI3
    Face3OrientE3 = -1
    Face3 = newnode(FaceList,
                    FaceTri_T(Face3Edge1, Face3Edge2, Face3Edge3, Face3OrientE1,
                              Face3OrientE2, Face3OrientE3))
    insertafter!(Face3, Face)

    Face.data.Edge1 = EdgeI1
    Face.data.OrientE1 = 1
    Face.data.Edge2 = EdgeI2
    Face.data.OrientE2 = 1
    Face.data.Edge3 = EdgeI3
    Face.data.OrientE3 = 1
end

function RefineFaceTriangularGrid!(TriangularGrid)
    EdgeList = TriangularGrid.EdgeList
    FaceList = TriangularGrid.FaceList

    Face = head(FaceList)
    while ~attail(Face)
        RefineFace!(Face, EdgeList, FaceList)
        Face = Face.next
        Face = Face.next
        Face = Face.next
        Face = Face.next
    end
end

function NumberingTriangularGrid!(TriangularGrid)
    NumNodes = 0
    NodeL = head(TriangularGrid.NodeList)
    while ~attail(NodeL)
        NumNodes += 1
        NodeL.data.Number = NumNodes
        NodeL = NodeL.next
    end
    TriangularGrid.NumNodes = NumNodes

    NumEdges = 0
    EdgeL = head(TriangularGrid.EdgeList)
    while ~attail(EdgeL)
        NumEdges += 1
        EdgeL.data.Number = NumEdges
        push!(EdgeL.data.Node1.data.Edge, NumEdges)
        push!(EdgeL.data.Node2.data.Edge, NumEdges)
        EdgeL = EdgeL.next
    end
    TriangularGrid.NumEdges = NumEdges

    NumFaces = 0
    FaceL = head(TriangularGrid.FaceList)
    while ~attail(FaceL)
        NumFaces += 1
        FaceL.data.Number = NumFaces
        push!(FaceL.data.Edge1.data.Face, NumFaces)
        push!(FaceL.data.Edge2.data.Face, NumFaces)
        push!(FaceL.data.Edge3.data.Face, NumFaces)
        FaceL = FaceL.next
    end
    TriangularGrid.NumFaces = NumFaces
end

function MidPoint(Face)
    P = Point()
    if Face.data.OrientE1 == 1
        P = Face.data.Edge1.data.Node1.data.P
    else
        P = Face.data.Edge1.data.Node2.data.P
    end
    if Face.data.OrientE2 == 1
        P = P + Face.data.Edge2.data.Node1.data.P
    else
        P = P + Face.data.Edge2.data.Node2.data.P
    end
    if Face.data.OrientE3 == 1
        P = P + Face.data.Edge3.data.Node1.data.P
    else
        P = P + Face.data.Edge3.data.Node2.data.P
    end
    M = Div(P, Norm(P))
end

function CircumCenterPoint(Face)
    P = Point()
    if Face.data.OrientE1 == 1
        P1 = Face.data.Edge1.data.Node1.data.P
    else
        P1 = Face.data.Edge1.data.Node2.data.P
    end
    if Face.data.OrientE2 == 1
        P2 = Face.data.Edge2.data.Node1.data.P
    else
        P2 = Face.data.Edge2.data.Node2.data.P
    end
    if Face.data.OrientE3 == 1
        P3 = Face.data.Edge3.data.Node1.data.P
    else
        P3 = Face.data.Edge3.data.Node2.data.P
    end
    M = CircumCenter(P1, P2, P3)
end

mutable struct Node
    P::Point
    N::Int
    NG::Int
    E::Array{Int, 1}
    F::Array{Int, 1}
    FG::Array{Int, 1}
    FP::Array{Int, 1}
    MasterSlave::Int
    Type::Char
end

function Point()
    return Point(0.0,
                 0.0,
                 0.0)
end

function Node()
    P = Point()
    N = 0
    NG = 0
    E = zeros(Int, 0)
    F = zeros(Int, 0)
    FG = zeros(Int, 0)
    FP = zeros(Int, 0)
    MasterSlave = 0
    Type = ' '
    return Node(P,
                N,
                NG,
                E,
                F,
                FG,
                FP,
                MasterSlave,
                Type)
end

function Node(Point::Point, Pos::Int, Type::Char)
    N = Node()
    N.P = Point
    N.N = Pos
    N.Type = Type
    return N
end

function FacesInNodes!(Nodes, Faces)
    NumNodes = size(Nodes, 1)
    NumFaces = size(Faces, 1)
    NumFacesPerNode = zeros(Int, NumNodes, 1)
    for iF in 1:NumFaces
        Face = Faces[iF]
        for iN in 1:size(Face.N, 1)
            NumFacesPerNode[Face.N[iN]] = NumFacesPerNode[Face.N[iN]] + 1
        end
    end
    FacesPerNode = zeros(Int, NumNodes, maximum(NumFacesPerNode))
    NumFacesPerNode = zeros(Int, NumNodes, 1)
    for iF in 1:NumFaces
        Face = Faces[iF]
        for iN in 1:size(Face.N, 1)
            NumFacesPerNode[Face.N[iN]] = NumFacesPerNode[Face.N[iN]] + 1
            FacesPerNode[Face.N[iN], NumFacesPerNode[Face.N[iN]]] = iF
        end
    end
    for iN in 1:NumNodes
        Nodes[iN].F = FacesPerNode[iN, 1:NumFacesPerNode[iN]]
    end
end

function AreaSphericalTriangle(P1, P2, P3)
    P1Loc = P1 / norm(P1)
    P2Loc = P2 / norm(P2)
    P3Loc = P3 / norm(P3)
    P1P2P3 = dot(P1Loc, P2Loc) + dot(P2Loc, P3Loc) + dot(P3Loc, P1Loc)
    P1_P2P3 = dot(P1Loc, cross(P2Loc, P3Loc))
    area = 2.0 * atan(abs(P1_P2P3) / (1.0 + P1P2P3))
end

function AreaFace(Face, Nodes)
    P1 = Nodes[Face.N[1]].P
    Area = 0.0
    for i in 2:(length(Face.N) - 1)
        P2 = Nodes[Face.N[i]].P
        P3 = Nodes[Face.N[i + 1]].P
        Area += AreaSphericalTriangle(P1, P2, P3)
    end
    return Area
end

function OrientFaceSphere(n::Point, m::Point)
    Orient = dot(n, m)
    return Orient
end

function TriangularGridToGrid(backend, FT, TriangularGrid, Rad, nz; ChangeOrient = 3)
    nBar = [0 1 1
            -1 1 0]
    Dim = 3
    Type = Tri()
    Form = SphericalGrid()

    NumNodes = TriangularGrid.NumNodes

    Nodes = map(1:NumNodes) do i
        Node()
    end

    NodeL = head(TriangularGrid.NodeList)
    NumNodes = 0
    while ~attail(NodeL)
        NumNodes += 1
        @show NodeL.data.P
        @show Rad
        Nodes[NumNodes] = Node(Rad * NodeL.data.P, NumNodes, 'N')
        NodeL = NodeL.next
    end

    NumEdges = TriangularGrid.NumEdges

    Edges = map(1:NumEdges) do i
        Edge()
    end

    EdgeL = head(TriangularGrid.EdgeList)
    NumEdges = 0
    while ~attail(EdgeL)
        NumEdges += 1
        n1 = EdgeL.data.Node1.data.Number
        n2 = EdgeL.data.Node2.data.Number
        Edges[NumEdges] = Edge(sort([n1; n2]), Nodes, NumEdges, NumEdges, "", NumEdges;
                               Form, Rad)
        EdgeL = EdgeL.next
    end

    NumFaces = TriangularGrid.NumFaces

    Faces = map(1:NumFaces) do i
        Face()
    end

    FaceL = head(TriangularGrid.FaceList)
    NumFaces = 0
    while ~attail(FaceL)
        NumFaces += 1
        e1 = FaceL.data.Edge1.data.Number
        e2 = FaceL.data.Edge2.data.Number
        e3 = FaceL.data.Edge3.data.Number
        s1 = sum(Edges[e1].N)
        s2 = sum(Edges[e2].N)
        s3 = sum(Edges[e3].N)
        permu = sortperm([s1; s2; s3])
        eee = [e1 e2 e3]
        ee = [eee[permu[1]]; eee[permu[3]]; eee[permu[2]]]
        (Faces[NumFaces], Edges) = Face(ee, Nodes, Edges, NumFaces, "Sphere",
                                        OrientFaceSphere;
                                        P = zeros(Float64, 0, 0), Form = Form, Rad = Rad,
                                        ChangeOrient = ChangeOrient)
        FaceL = FaceL.next
    end
    NumNodes = size(Nodes, 1)
    NumEdges = size(Edges, 1)
    NumFaces = size(Faces, 1)
    NumEdgesI = size(Edges, 1)

    FacesInNodes!(Nodes, Faces)

    zP = zeros(nz)
    z = zeros(FT, nz + 1)
    dzeta = zeros(nz)
    H = 0.0
    NumFacesB = 0
    NumFacesG = 0
    NumEdgesB = 0
    NumEdgesG = 0
    NumNodesB = 0
    NumNodesG = 0
    nBar3 = zeros(0, 0)
    AdaptGrid = ""
    EF = zeros(Int, 0, 0)
    FE = zeros(Int, 0, 0)

    return GridStruct{FT,
                      typeof(EF),
                      typeof(z)}(nz,
                                 zP,
                                 z,
                                 dzeta,
                                 H,
                                 NumFaces,
                                 NumFacesB,
                                 NumFacesG,
                                 Faces,
                                 NumEdges,
                                 NumEdgesB,
                                 NumEdgesG,
                                 Edges,
                                 NumNodes,
                                 NumNodesB,
                                 NumNodesG,
                                 Nodes,
                                 Form,
                                 Type,
                                 Dim,
                                 Rad,
                                 nBar3,
                                 nBar,
                                 AdaptGrid,
                                 EF,
                                 FE)
end

mutable struct Face
    N::Array{Int, 1}
    E::Array{Int, 1}
    F::Int
    FG::Int
    n::Point
    Mid::Point
    OrientE::Array{Int, 1}
    Type::String
    P::Array{Point, 1}
    Stencil::Array{Int, 1}
    Area::Float64
    Radius::Float64
    Orientation::Int64
end

function Face()
    N = zeros(Int, 0)
    E = zeros(Int, 0)
    F = 0
    FG = 0
    n = Point()
    Mid = Point()
    OrientE = zeros(Int, 0)
    Type = ""
    P = Array{Point}(undef, 0)
    Stencil = zeros(Int, 0)
    Area::Float64 = 0
    Radius::Float64 = 0
    Orientation::Int = 1
    return Face(N,
                E,
                F,
                FG,
                n,
                Mid,
                OrientE,
                Type,
                P,
                Stencil,
                Area,
                Radius,
                Orientation)
end

function Face(EdgesF::Array{Int, 1}, Nodes, Edges, Pos, Type, OrientFace;
              Form = CartesianGrid(), Rad = 1.0,
              P::Array{Float64, 2} = [], ChangeOrient = 3, MidFace = nothing)
    F = Face()
    if EdgesF[1] == 0
        return (F, Edges)
    end

    nE = size(EdgesF, 1)
    F.F = Pos
    F.Type = Type
    # TODO: check translation
    @inbounds for iE in 1:nE
        Edges[EdgesF[iE]].NumF += 1
        Edges[EdgesF[iE]].F[Edges[EdgesF[iE]].NumF] = Pos
    end
    #Sort edges
    F.E = zeros(Int, nE)
    F.E[1] = EdgesF[1]
    N2 = Edges[F.E[1]].N[2]
    @inbounds for iE in 2:nE
        @inbounds for iE1 in iE:nE
            if N2 == Edges[EdgesF[iE1]].N[1]
                F.E[iE] = EdgesF[iE1]
                N2 = Edges[EdgesF[iE1]].N[2]
                EdgesF[iE1] = EdgesF[iE]
                EdgesF[iE] = F.E[iE]
                break
            elseif N2 == Edges[EdgesF[iE1]].N[2]
                F.E[iE] = EdgesF[iE1]
                N2 = Edges[EdgesF[iE1]].N[1]
                EdgesF[iE1] = EdgesF[iE]
                EdgesF[iE] = F.E[iE]
                break
            end
        end
    end
    F.N = zeros(Int, nE)
    F.N[1:2] = Edges[F.E[1]].N
    @inbounds for iE in 2:(nE - 1)
        if F.N[iE] == Edges[F.E[iE]].N[1]
            F.N[iE + 1] = Edges[F.E[iE]].N[2]
        else
            F.N[iE + 1] = Edges[F.E[iE]].N[1]
        end
    end
    if P == zeros(Float64, 0, 0)
        F.P = Array{Point}(undef, size(F.N, 1))
        @inbounds for i in 1:size(F.N, 1)
            F.P[i] = Nodes[F.N[i]].P
        end
    else
        F.P = Array{Point}(undef, size(F.N, 1))
        @inbounds for i in 1:size(F.N, 1)
            F.P[i] = Point(P[:, i])
        end
    end
    if Form == SphericalGrid()
        F.Area = AreaFace(F, Nodes) * Rad * Rad
    else
        PT = Point([0.0, 0.0, 0.0])
        @inbounds for i in 1:(nE - 1)
            PT = PT + cross(F.P[i], F.P[i + 1])
        end
        PT = PT + cross(F.P[nE], F.P[1])
        F.Area = 0.5 * norm(PT)
    end
    if MidFace === nothing
        @inbounds for i in 1:nE
            F.Mid = F.Mid + F.P[i]
        end
        F.Mid = F.Mid / Float64(nE)
    else
        F.Mid.x = MidFace.x
        F.Mid.y = MidFace.y
        F.Mid.z = MidFace.z
    end
    if Form == SphericalGrid()
        F.Mid = F.Mid / norm(F.Mid) * Rad
    end
    if Form == SphericalGrid()
        F.Radius = Rad
    else
    end

    NumE = size(EdgesF, 1)
    F.n = cross(F.P[NumE], F.P[1])
    @inbounds for i in 1:(NumE - 1)
        F.n = F.n + cross(F.P[i], F.P[i + 1])
    end
    F.n = F.n / norm(F.n)
    if OrientFace(F.n, F.Mid) < 0
        F.Orientation = -1
        if NumE > ChangeOrient
            #Change Orientation
            NTemp = copy(F.N)
            ETemp = copy(F.E)
            PTemp = copy(F.P)
            @inbounds for i in 1:nE
                F.N[i] = NTemp[nE - i + 1]
                F.P[i] = PTemp[nE - i + 1]
            end
            @inbounds for i in 1:(nE - 1)
                F.E[i] = ETemp[nE - i]
            end
            F.n = -F.n
            F.Orientation = 1
        end
    end
    F.OrientE = zeros(Int, NumE)
    for i in 1:NumE
        iE = F.E[i]
        if Edges[iE].N[1] == F.N[i]
            F.OrientE[i] = F.Orientation
            F.OrientE[i] = 1
        else
            F.OrientE[i] = -F.Orientation
            F.OrientE[i] = -1
        end
    end
    return F, Edges
end

function DelaunayGridToPolyGrid(backend, FT, TriangularGrid, Rad, nz)
    nBar = [0 1 0 1
            -1 0 -1 0]
    Dim = 3
    Type = Tri()
    Rad = Rad
    Form = SphericalGrid()

    NumNodes = TriangularGrid.NumFaces

    Nodes = map(1:NumNodes) do i
        Node()
    end

    FaceL = head(TriangularGrid.FaceList)
    NumNodes = 0
    while ~attail(FaceL)
        NumNodes += 1
        PM = MidPoint(FaceL)
        PC = CircumCenterPoint(FaceL)
        if dot(PM, PC) < 0
            PC = -PC
        end
        Nodes[NumNodes] = Node(Rad * PC, NumNodes, 'N')
        FaceL = FaceL.next
    end

    NumEdges = TriangularGrid.NumEdges

    Edges = map(1:NumEdges) do i
        Edge()
    end

    EdgeL = head(TriangularGrid.EdgeList)
    NumEdges = 0
    while ~attail(EdgeL)
        NumEdges += 1
        n1 = EdgeL.data.Face[1]
        n2 = EdgeL.data.Face[2]
        Edges[NumEdges] = Edge([n1, n2], Nodes, NumEdges, NumEdges, "", NumEdges; Form, Rad)
        EdgeL = EdgeL.next
    end

    NumFaces = TriangularGrid.NumNodes

    Faces = map(1:NumFaces) do i
        Face()
    end

    NodeL = head(TriangularGrid.NodeList)
    NumFaces = 0
    while ~attail(NodeL)
        NumFaces += 1
        e = NodeL.data.Edge
        (Faces[NumFaces], Edges) = Face(e, Nodes, Edges, NumFaces, "Sphere",
                                        OrientFaceSphere;
                                        Form = Form, Rad = Rad, P = zeros(Float64, 0, 0))
        NodeL = NodeL.next
    end
    NumNodes = size(Nodes, 1)
    NumEdges = size(Edges, 1)
    NumFaces = size(Faces, 1)
    NumEdgesI = size(Edges, 1)

    NumEdgesB = 0

    FacesInNodes!(Nodes, Faces)
    SortFacesInNodes!(Nodes, Faces)

    TestOrientation(Faces, Nodes)

    zP = zeros(nz)
    z = KernelAbstractions.zeros(backend, FT, nz + 1)
    dzeta = zeros(nz)
    H = 0.0
    NumFacesB = 0
    NumFacesG = 0
    NumEdgesB = 0
    NumEdgesG = 0
    NumNodesB = 0
    NumNodesG = 0
    nBar3 = zeros(0, 0)
    nBar = zeros(0, 0)
    AdaptGrid = ""
    EF = KernelAbstractions.zeros(backend, Int, 0, 0)
    FE = KernelAbstractions.zeros(backend, Int, 0, 0)

    return GridStruct{FT,
                      typeof(EF),
                      typeof(z)}(nz,
                                 zP,
                                 z,
                                 dzeta,
                                 H,
                                 NumFaces,
                                 NumFacesB,
                                 NumFacesG,
                                 Faces,
                                 NumEdges,
                                 NumEdgesB,
                                 NumEdgesG,
                                 Edges,
                                 NumNodes,
                                 NumNodesB,
                                 NumNodesB,
                                 Nodes,
                                 Form,
                                 Type,
                                 Dim,
                                 Rad,
                                 nBar3,
                                 nBar,
                                 AdaptGrid,
                                 EF,
                                 FE)
end

function TriangularGrid(backend, FT, RefineLevel, RadEarth, nz; ChangeOrient = 3)
    IcosahedronGrid = CreateIcosahedronGrid()
    for iRef in 1:RefineLevel
        RefineEdgeTriangularGrid!(IcosahedronGrid)
        RefineFaceTriangularGrid!(IcosahedronGrid)
    end
    NumberingTriangularGrid!(IcosahedronGrid)
    Grid = TriangularGridToGrid(backend, FT, IcosahedronGrid, RadEarth, nz;
                                ChangeOrient = ChangeOrient)
    # OrientTriangle(Grid)
    return Grid
end

function DelaunayGrid(backend, FT, RefineLevel, RadEarth, nz)
    IcosahedronGrid = Grids.CreateIcosahedronGrid()
    for iRef in 1:RefineLevel
        Grids.RefineEdgeTriangularGrid!(IcosahedronGrid)
        Grids.RefineFaceTriangularGrid!(IcosahedronGrid)
    end
    Grids.NumberingTriangularGrid!(IcosahedronGrid)
    Grids.DelaunayGridToPolyGrid(backend, FT, IcosahedronGrid, RadEarth, nz)
end

mutable struct Edge
    N::Array{Int, 1}
    E::Int
    EG::Int
    EI::Int
    ET::Int
    NumF::Int
    F::Array{Int, 1}
    FG::Array{Int, 1}
    FP::Array{Int, 1}
    t::Point
    n::Point
    a::Float64
    Mid::Point
    FE::Array{Int, 1}
    Type::String
    MasterSlave::Int
end

"""
  Edge()

This is my documentation
"""
function Edge()
    N = zeros(Int, 2)
    E = 0
    EG = 0
    EI = 0
    ET = 0
    NumF = 0
    F = zeros(Int, 2)
    FG = zeros(Int, 2)
    FP = zeros(Int, 2)
    t = Point()
    n = Point()
    a = 0.0
    Mid = Point()
    FE = zeros(Int, 0)
    Type = ""
    MasterSlave = 0
    return Edge(N,
                E,
                EG,
                EI,
                ET,
                NumF,
                F,
                FG,
                FP,
                t,
                n,
                a,
                Mid,
                FE,
                Type,
                MasterSlave)
end

mutable struct GreatCircle
    P1::Point
    P2::Point
end

function SizeGreatCircle(C::GreatCircle)
    return acos(dot(C.P1, C.P2) / (norm(C.P1) * norm(C.P2)))
end

function SizeGreatCircle(Lon1, Lat1, Lon2, Lat2)
    return acos(sin(Lat1) * sin(Lat2) +
                cos(Lat1) * cos(Lat2) * cos(Lon2 - Lon1))
end

function SizeGreatCircle(Edge::Edge, Nodes)
    P1 = Nodes[Edge.N[1]].P
    P2 = Nodes[Edge.N[2]].P
    return acos(dot(P1, P2) / (norm(P1) * norm(P2)))
end

function SizeGreatCircle(P1::Point, P2::Point)
    return acos(dot(P1, P2) / (norm(P1) * norm(P2)))
end

function Edge(NodesE, Nodes, PosG, PosI, Type, PosT = nothing; Form = CartesianGrid(),
              Rad = 1.0)
    E = Edge()
    E.E = PosG
    E.EI = PosI
    if PosT ≠ nothing
        E.ET = PosT
    else
        E.ET = 0
    end
    E.N = NodesE
    E.t = Nodes[E.N[2]].P - Nodes[E.N[1]].P
    if Form == SphericalGrid()
        E.a = SizeGreatCircle(Nodes[E.N[2]].P, Nodes[E.N[1]].P) * Rad
    else
        E.a = norm(E.t)
    end
    E.t = E.t / norm(E.t)
    E.Mid = 0.5 * (Nodes[E.N[1]].P + Nodes[E.N[2]].P)
    if Form == SphericalGrid()
        E.Mid = E.Mid * (Rad / norm(E.Mid))
    end
    k = zeros(3)
    t = zeros(3)
    n = zeros(3)
    k[1] = E.Mid.x
    k[2] = E.Mid.y
    k[3] = E.Mid.z
    k = k / norm(k)
    t[1] = E.t.x
    t[2] = E.t.y
    t[3] = E.t.z
    n = LinearAlgebra.cross(t, k)
    E.n.x = n[1]
    E.n.y = n[2]
    E.n.z = n[3]
    E.Type = Type
    return E
end

mutable struct GridStruct{FT <: AbstractFloat,
                          IT2 <: AbstractArray,
                          AT1 <: AbstractArray}
    nz::Int
    zP::Array{FT, 1}
    z::AT1
    dzeta::Array{FT, 1}
    H::FT
    NumFaces::Int
    NumFacesB::Int
    NumFacesG::Int
    Faces::Array{Face, 1}
    NumEdges::Int
    NumEdgesB::Int
    NumEdgesG::Int
    Edges::Array{Edge, 1}
    NumNodes::Int
    NumNodesB::Int
    NumNodesG::Int
    Nodes::Array{Node, 1}
    Form::GridForm
    Type::ElementType
    Dim::Int
    Rad::FT
    nBar3::Array{FT, 2}
    nBar::Array{FT, 2}
    AdaptGrid::Any
    EF::IT2
    FE::IT2
end

function AddVerticalGrid!(Grid::GridStruct, nz::Int, H)
    Grid.zP = zeros(nz)
    z = zeros(nz + 1)
    Grid.dzeta = zeros(nz)
    Grid.H = H
    @. Grid.dzeta = H / nz
    Grid.H = H
    for i in 2:(nz + 1)
        z[i] = z[i - 1] + Grid.dzeta[i - 1]
    end
    for i in 1:nz
        Grid.dzeta[i] = z[i + 1] - z[i]
        Grid.zP[i] = 0.5 * (z[i] + z[i + 1])
    end
    copyto!(Grid.z, z)
end

mutable struct MetricDGStruct{FT <: AbstractFloat,
                              AT2 <: AbstractArray,
                              AT3 <: AbstractArray,
                              AT4 <: AbstractArray,
                              AT5 <: AbstractArray,
                              AT6 <: AbstractArray}
    J::AT4
    X::AT5
    dXdxI::AT6
    Rotate::AT6
    dz::AT2
    zP::AT2
    xS::AT2
    VolSurfH::AT4
    NH::AT5
    VolSurfV::AT3
    NV::AT4
end

function MetricCreate(backend, FT, nQuad, OPZ, NF, nz, NumG, placeholder)
    J = KernelAbstractions.zeros(backend, FT, nQuad, OPZ, nz, NF)
    X = KernelAbstractions.zeros(backend, FT, nQuad, OPZ, 3, nz, NF)
    dXdxI = KernelAbstractions.zeros(backend, FT, 3, 3, OPZ, nQuad, nz, NF)
    Rotate = KernelAbstractions.zeros(backend, FT, 3, 3, OPZ, nQuad, nz, NF)
    dz = KernelAbstractions.zeros(backend, FT, 0, 0)
    zP = KernelAbstractions.zeros(backend, FT, 0, 0)
    xS = KernelAbstractions.zeros(backend, FT, 2, NumG)
    VolSurfH = KernelAbstractions.zeros(backend, FT, 0, 0, 0, 0)
    NH = KernelAbstractions.zeros(backend, FT, 0, 0, 0, 0, 0)
    VolSurfV = KernelAbstractions.zeros(backend, FT, 0, 0, 0)
    NV = KernelAbstractions.zeros(backend, FT, 0, 0, 0, 0)
    return MetricDGStruct{FT,
                          typeof(zP),
                          typeof(VolSurfV),
                          typeof(J),
                          typeof(X),
                          typeof(dXdxI)}(J,
                                         X,
                                         dXdxI,
                                         Rotate,
                                         dz,
                                         zP,
                                         xS,
                                         VolSurfH,
                                         NH,
                                         VolSurfV,
                                         NV)
end

mutable struct MetricCGStruct{FT <: AbstractFloat,
                              AT2 <: AbstractArray,
                              AT3 <: AbstractArray,
                              AT4 <: AbstractArray,
                              AT5 <: AbstractArray,
                              AT6 <: AbstractArray}
    J::AT4
    X::AT5
    dXdxI::AT6
    dXdx::AT6
    nSS::AT2
    nS::AT3
    FS::AT2
    dz::AT2
    zP::AT2
    JC::AT3
    JCW::AT3
    xS::AT2
    M::AT3
    MMass::AT2
end

function MetricCompute(FT, FE, Grid, NumberThreadGPU, zS)
    DoF = FE.DoF
    DoFE = FE.DoFE
    M = FE.OrdPolyZ + 1
    NF = Grid.NumFaces
    NE = Grid.NumEdges
    Nz = Grid.nz

    EdgeFace!(Grid)
    Metric = MetricCreate(FT, DoF, M, NF, Nz, FE.NumG, FE)
    Metric.zP = zeros(FT, Nz, FE.NumG)
    Metric.dz = zeros(FT, Nz, FE.NumG)
    F = PointsFromGrid(FT, Grid)
    FillX!(Metric, FE, F, Grid, zS)
    FillRotate!(Metric, FE, Grid)
    FillContravariant!(Metric, FE, Grid, Grid.Type, Model.MetricType)
    FillDet!(Metric, FE, Grid)
    GridSizeDGKernel!(FE, Metric, Grid.Rad, NumberThreadGPU, Grid.Form)
    MetricLowerBoundary!(Metric, FE, Grid, NumberThreadGPU, Grid.Form)
    NormalH!(Metric, FE, Grid, NumberThreadGPU, Grid.Type)
    NormalV!(Metric, FE, Grid, NumberThreadGPU)

    return Metric
end

ChangeOrient = 3
RadEarth = EARTH_RADIUS
nz = 4
RefineLevel = 1
FT = Float64
Grid = TriangularGrid(nothing, FT, RefineLevel, RadEarth, nz; ChangeOrient = ChangeOrient)
H = 30000
AddVerticalGrid!(Grid, nz, H)
NumberThreadGPU = 4
Metric = MetricCompute(FT, Grid.FE, Grid, NumberThreadGPU, zS)
