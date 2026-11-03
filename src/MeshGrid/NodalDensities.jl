"""
Functions for calculating nodal densities (mean or least squares).

This part provides functions to transform element-based values to nodal values
using different approaches: simple averaging or least squares interpolation.
"""

# Compute nodal values by averaging element values connected to each node
function elementToNodalValues(
    mesh::Mesh,
    elVals::Vector{Float64},
)::Vector{Float64}
    nodalVals = zeros(Float64, mesh.nnp)

    for A = 1:mesh.nnp
        count = 0
        for el in mesh.INE[A]
            nodalVals[A] += elVals[el]
            count += 1
        end
        nodalVals[A] /= count
    end

    return nodalVals
end
##################################################################

## Store nodal coordinates for each element
mutable struct NodalCoordinatesInElement{T<:Array}
    x::T
    y::T
    z::T
end

"""
    NodePosition3D(mesh::AbstractMesh)

Extract node coordinates for each element.
Returns a NodalCoordinatesInElement struct containing 3D coordinate arrays
where each row corresponds to an element ID.
"""
function NodePosition3D(mesh::AbstractMesh)
    table = size(mesh.IEN)
    EN_x = zeros(Float64, table)
    EN_y = zeros(Float64, table)
    EN_z = zeros(Float64, table)

    for i = 1:(length(mesh.IEN[:, 1]))
        for j = 1:(length(mesh.IEN[1, :]))
            ID_n = Int(mesh.IEN[i, j])
            x = mesh.X[1, ID_n]
            y = mesh.X[2, ID_n]
            z = mesh.X[3, ID_n]
            EN_x[i, j] = x
            EN_y[i, j] = y
            EN_z[i, j] = z
        end
    end
    EN = NodalCoordinatesInElement(EN_x, EN_y, EN_z)

    return EN # Output consists of three arrays - coordinates of nodes for each element (row position = element ID)
end


"""
    GeometricCentre(mesh::AbstractMesh, EN::NodalCoordinatesInElement)

Calculate geometric centers of elements (equivalent to Gauss points for first order elements).
Returns a matrix where each row contains the 3D coordinates of an element's center.
"""
function GeometricCentre(mesh::AbstractMesh, EN::NodalCoordinatesInElement)
    Centre = zeros(mesh.nel, length(mesh.X[:, 1]))
    for i = 1:mesh.nel
        G_x = mean(EN.x[:, i]) # G_x = dot(EN_x[:, i], N)
        G_y = mean(EN.y[:, i]) # G_y = dot(EN_y[:, i], N)
        G_z = mean(EN.z[:, i]) # G_z = dot(EN_z[:, i], N)
        Centre[i, :] = [G_x, G_y, G_z]
    end
    return Centre
end


"""
    DenseInNodes(mesh::Mesh, ρ::Vector)

Calculate density at nodes based on element densities using different
interpolation methods depending on how many elements share each node.
"""
function DenseInNodes(
    mesh::Mesh, 
    ρ::Vector,
)
    EN = NodePosition3D(mesh)
    Centre = GeometricCentre(mesh, EN)

    dense_nodes = zeros(Float64, mesh.nnp) # Number of nodes
    for i = 1:mesh.nnp # Loop through nodes (from ID_1)
        No1E = length(mesh.INE[i])
        if No1E == 1 # Node belongs to one element only (corners)
            dense_nodes[i] = ρ[Int((mesh.INE[i])[1])]
        elseif No1E > 1 && No1E < 4
            dense_nodes[i] = FilterForNodalDensity(i, Centre, mesh, ρ)
        elseif No1E > 3 # More elements (4+) -> least squares fitting (weighted average on the boundary)
            dense_nodes[i] = NodalDensityLeastSquares(i, Centre, mesh, ρ)
        end
    end
    return dense_nodes
end


"""
    FilterForNodalDensity(i::Int, Centre::Matrix, mesh::Mesh, ρ::Vector)

Calculate nodal density using a distance-weighted approach for nodes
that belong to 2 or 3 elements. Also used for boundary nodes with 4+ elements,
where the least squares fit is not reliable.

The nodal density is the weighted average of the element densities with the weights
1 - Lⱼ / (1.2 L_max), where Lⱼ is the distance from the node to the centre of element j and
L_max the largest of them. The weights fall linearly with the distance (they are not
inverse-distance weights), the factor 1.2 keeps the farthest element at the weight 1/6.
"""
function FilterForNodalDensity(
    i::Int,
    Centre::Matrix,
    mesh::Mesh,
    ρ::Vector,
)
    L = zeros(length(mesh.INE[i]))

    for j = 1:length(mesh.INE[i])
        L[j] = norm(mesh.X[:, i] - Centre[mesh.INE[i][j], :])
    end
    Lmax = maximum(L) * 1.2
    dm = 0.0
    deleni = 0.0
    for j = 1:length(mesh.INE[i])
        dm += ρ[Int((mesh.INE[i])[j])] * (1 - L[j] / Lmax)
        deleni += (1 - L[j] / Lmax)
    end
    return DN = dm / deleni
end


"""
    NodalDensityLeastSquares(i::Int, Centre::Matrix, mesh::Mesh, ρ::Vector)

Calculate nodal density using least squares fitting for nodes that belong to 4 or more elements.
A linear function is fitted through the densities in the element centres. The coordinates
are measured from the node, so the nodal density is simply the first coefficient a₁.
Nodes on the boundary of the domain use the weighted average (FilterForNodalDensity) instead.
"""
function NodalDensityLeastSquares(
    i::Int,
    Centre::Matrix,
    mesh::Mesh,
    ρ::Vector,
)
    nnp_i = length(mesh.INE[i])
    xₙ = mesh.X[:, i]

    εₘ₁ = 1.e-8 # limit for "centres lie in one plane"
    εₘ₂ = 2.0   # limit for "node is too far from the middle of the centres"

    # Mean distance from the node to the element centres
    L_sum = 0.0
    for j = 1:nnp_i
        L_sum += norm(Centre[mesh.INE[i][j], :] - xₙ)
    end
    L = L_sum / nnp_i

    A = zeros(nnp_i, mesh.nsd + 1)
    b = zeros(nnp_i)
    for j = 1:nnp_i
        el = mesh.INE[i][j]
        A[j, :] = [1; (Centre[el, :] - xₙ) / L] # centre position relative to the node
        b[j] = ρ[el]
    end

    # Singular value decomposition A = U * Diagonal(S) * Vᵀ
    # S (singular values) show how well the centres span 3D, F \ b gives the least squares fit
    F = svd(A)

    # Centres lie in one plane (node on a face of the domain) -> fit is not unique
    if minimum(F.S) <= εₘ₁ * maximum(F.S)
        return FilterForNodalDensity(i, Centre, mesh, ρ)
    end

    # Leverage: 1 = node in the middle of the centres, larger = node off to the side
    leverage = nnp_i * sum((F.V[1, :] ./ F.S) .^ 2)
    if leverage > εₘ₂ # fit would extrapolate -> use weighted average
        return FilterForNodalDensity(i, Centre, mesh, ρ)
    end

    a = F \ b # least squares fit
    DN = a[1] # fit value at the node
    return DN
end

"""
    KeepRange(DN::Float64)

Constrain density value to stay within valid range [0, 1].
"""
function KeepRange(DN::Float64)
    if DN < 0
        DN = 0.0
    elseif DN > 1
        DN = 1
    end
    return DN
end
