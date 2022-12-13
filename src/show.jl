#---Helper functions used to print common features (e.g. basis vectors)---------------------------#
"""
    Xtal.vector_string(v::AbstractVector{<:Real}; brackets=true) -> String

Prints a representation of a vector as a string.
"""
function vector_string(v::AbstractVector{<:Real}; brackets=true)
    # Format the numbers within a vector
    tostr(x, n) = lpad(@sprintf("%f", x), n)
    return "["^brackets * join(tostr.(v, 10)) *" ]"^brackets
end

"""
    Xtal.basis_string(M::AbstractMatrix{<:Real}; letters=true) -> Vector{String}

Prints each basis vector with an associated letter.

# Example
```
a: [  0.000000  3.500000  3.500000 ]   (4.949747 Å)
b: [  3.500000  0.000000  3.500000 ]   (4.949747 Å)
c: [  3.500000  3.500000  0.000000 ]   (4.949747 Å)
```
"""
function basis_string(
    M::AbstractMatrix{<:Real};
    pad="  ",
    brackets=true,
    letters=true,
    length=true,
    unit=""
)
    # Format the numbers within a vector
    tostr(x, n) = lpad(@sprintf("%f", x), n)
    # Letters are generated by incrementing up from character 0x60
    # Char(0x61) = 'a', Char(0x62) = 'b'...
    # Letters should work up to 26 dimensions, but who's gonna deal with 26D crystals?
    # Bosonic string theorists, maybe?
    return [
        pad * string(Char(0x60 + n), ':', ' ')^letters *
        vector_string(M[:,n], brackets=brackets) *
        ("   (" * tostr(norm(M[:,n]), 0))^length * " "^!isempty(unit) * unit * ")"
        for n in 1:size(M)[2]
    ]
end

basis_string(b::RealBasis, kwargs...) = basis_string(matrix(b), unit="Å", kwargs...)
basis_string(b::ReciprocalBasis, kwargs...) = basis_string(matrix(b), unit="Å⁻¹", kwargs...)

"""
    Xtal.printbasis(io::IO, b; kwargs...)

Prints the result of `basis_string()` to `io`.
"""
function printbasis(io::IO, M::AbstractMatrix{<:Real}; letters=true, unit="", pad=0)
    s = basis_string(M, letters=letters, unit=unit)
    print(io, join(" "^pad .* s, "\n"))
end

printbasis(io::IO, b::RealBasis; kwargs...) = printbasis(io, matrix(b), unit="Å"; kwargs...)
printbasis(io::IO, b::ReciprocalBasis; kw...) = printbasis(io, matrix(b), unit="Å⁻¹"; kw...)
printbasis(io::IO, a; kwargs...) = printbasis(io, basis(a); kwargs...)

"""
    atom_string(a::AtomPosition; name=true, num=true)

Generates a string describing an atom position.
"""
function atom_string(a::AtomPosition; name=true, num=true, entrysz=4)
    # Format the numbers within a vector
    tostr(x) = lpad(@sprintf("%f", x), 10)
    return rpad(string(a.num), entrysz)^num * rpad(a.name, entrysz)^name * vector_string(a.pos)
end

# TODO: make this work with occupancy information, since we'll probably need that
function formula_string(v::AbstractVector{<:AtomPosition}, reduce=true)
    # Number of each type of atoms is stored in this vector
    atomcount = zeros(Int, size(ELEMENTS)...)
    # Get all atomic numbers in the 
    atomnos = [a.num for a in v]
    # Get all of the types of atoms
    for atom in atomnos
        # Skip dummy atoms
        atom == 0 && continue
        # Increment the atom counts
        atomcount[atom] += 1
    end
    # Create output string
    str = ""
    # Loop through all the atom counts
    for (atom, ct) in enumerate(atomcounts)
        # Skip any zeros
        ct == 0 && continue
        str *= ELEMENT_LOOKUP[atom] * string(ct) * space
    end
    return str
end

"""
    formula_string(a::AtomList; reduce=true) -> String

Generates a string giving the atomic formula for an `AtomicPosition`. By default, common factors
will be reduced.
"""
formula_string(a::AtomList; reduce=true) = formula_string(a.coord, reduce=reduce)

#---Actual show methods---------------------------------------------------------------------------#

# These are what you see when something is returned in the REPL.
# To define these methods for a type, just overload show(::IO, ::MIME"text/plain", ::T)
# To get the result as a string, just use repr("text/plain", x)

#---Types from lattices.jl (RealBasis, ReciprocalBasis)-------------------------------------------#

function Base.show(io::IO, ::MIME"text/plain", b::AbstractBasis)
    println(io, typeof(b), ":")
    printbasis(io, b, pad=2)
end

#---Types from atoms.jl (AtomPosition, AtomList)--------------------------------------------------#

function Base.show(io::IO, ::MIME"text/plain", a::AtomPosition; kwargs...)
    println(io, typeof(a), ":")
    println(io, "  ", atom_string(a; kwargs...))
end

function Base.show(io::IO, ::MIME"text/plain", a::AtomList; letters=true, kwargs...)
    # Print type name
    println(io, typeof(a), ":")
    # Print atomic positions
    println(io, "  ", natom(a), " atomic positions:")
    for atom in a.coord
        println(io, "    ", atom_string(atom; kwargs...))
    end
    # Print basis vectors
    println("  defined in terms of basis vectors:")
    printbasis(io, a, pad=2)
end

#---Types from data/realspace.jl------------------------------------------------------------------#

function Base.show(io::IO, ::MIME"text/plain", g::RealSpaceDataGrid)
    dimstring = join(string.(size(g)), "×") * " "
    println(io, dimstring, typeof(g), " with real space basis vectors:")
    print(io, join(basis_string(basis(g)), "\n"))
    @printf(io, "\nCell volume: %16.10f Å", volume(g))
    @printf(io, "\nVoxel size:  %16.10f Å", voxelsize(g))
end

#---Types from data/reciprocalspace.jl------------------------------------------------------------#

function Base.show(io::IO, ::MIME"text/plain", g::HKLData)
    dimstring = join(string.(size(g)), "×") * " "
    println(io, dimstring, typeof(g), " with reciprocal space basis vectors:")
    print(io, join(basis_string(basis(g), unit="Å⁻¹"), "\n"))
    print(io, "\nMaximum spatial frequencies:")
    for n in 1:length(basis(g))
        print(io, "\n  ", '`' + n, ":", )
        @printf(io, "%12.6f Å⁻¹", size(g)[n] * lengths(basis(g))[n])
    end
end

function Base.show(io::IO, ::MIME"text/plain", wf::ReciprocalWavefunction)
    println(io,
        typeof(wf), " with ",
        string(nspin(wf)), " spin", "s"^(nspin(wf) != 1), ", ",
        string(nkpt(wf)), " k-point", "s"^(nkpt(wf) != 1), ", and ",
        string(nband(wf)), " band", "s"^(nband(wf) != 1)
    )
    println(io, "Reciprocal space basis vectors:")
    print(io, join(basis_string(basis(wf)), "\n"))
end

#---Types from data/atomic.jl---------------------------------------------------------------------#

function Base.show(
    io::IO,
    ::MIME"text/plain",
    s::SphericalHarmonic{Lmax};
    showto = 3
) where Lmax
    # Don't include generated second type parameter
    print(io, "SphericalHarmonic{$Lmax}", ":\n", " "^13)
    # Only print up to l=3 components by default (kw showto)
    Lmax_eff = min(showto, Lmax)
    for m = -Lmax_eff:Lmax_eff
        print(io, rpad("m = $m", 12))
    end
    for l in 0:Lmax_eff
        print(io, "\n", rpad("l = $l:", 8))
        for m = -Lmax_eff:Lmax_eff
            if abs(m) <= l
                print(io, lpad(@sprintf("%6f", s[l,m]), 12))
            else
                print(io, " "^12)
            end
        end
    end
    if showto < Lmax
        print(io, "\n(higher order components omitted for brevity)")
    end
end

#---Types from crystals.jl------------------------------------------------------------------------#

function Base.show(io::IO, ::MIME"text/plain", xtal::Crystal{D}) where D
    println(io, typeof(xtal), " (space group ", xtal.sgno, "): ")
    # Print basis vectors
    println(io, "\n  Primitive basis vectors:")
    printbasis(io, xtal, pad=2, unit="Å")
    if xtal.transform != SMatrix{D,D,Float64}(LinearAlgebra.I)
        println(io, "\n\n  Conventional basis vectors:")
        printbasis(io, basis(xtal) * xtal.transform, pad=2, unit="Å")
    end
    # TODO: Add in more info about atomic positions, space group
    println(io, "\n\n  ", length(xtal.atoms), " atomic positions:")
    println(io, "    Num   ", "Name  ", "Position")
    for atom in xtal.atoms
        println(io, "    ", atom_string(atom, name=true, num=true, entrysz=6))
    end
    # Determine what basis the atomic coordinates are given in.
    # If the basis is zero, assume Cartesian coordinates in Å
    if basis(xtal.atoms) == zeros(RealBasis{D})
        print("  in Cartesian coordinates (assumed to be Å)")
    else
        print("  in fractional coordinates")
    end
end

function Base.show(io::IO, ::MIME"text/plain", x::CrystalWithDatasets)
    println(io, typeof(x), " containing:\n")
    show(io, MIME("text/plain"), x.xtal)
    print("\n\nand a ")
    show(io, MIME("text/plain"), x.data)
end

#---Other internal types--------------------------------------------------------------------------#

function Base.show(io::IO, ::MIME"text/plain", h::ABINITHeader)  
    println(io, "abinit ", repr(h.codvsn)[3:end-1], " header (version ", h.headform, "):")
    for name in fieldnames(ABINITHeader)[3:end]
        print(io, "  ", rpad(string(name) * ":", 18))
        x = getfield(h, name)
        if typeof(x) <: Union{<:Number,<:SVector}
            println(io, x)
        elseif typeof(x) <: SMatrix
            println(io, replace(repr("text/plain", x), "\n " => "\n" * " "^21))
        elseif typeof(x) <: AbstractArray
            sizestr = if length(size(x)) == 1
                lpad(string(length(x), "-element "), 12)
            else    
                join(string.(size(x)), "×") * " "
            end
            println(io, sizestr, typeof(x), "...")
        else
            println(io, typeof(x), "...")
        end
    end
end
