module CliffordAlgebrasPrettyTablesExt

using CliffordAlgebras
import PrettyTables: pretty_table

# Provide a method for the internal renderer that delegates to PrettyTables
function CliffordAlgebras._render_table(
    io::IO,
    table;
    show_row_number::Bool = false,
    crop = :none,
    show_header::Bool = false,
    alignment = :c,
    vlines = nothing,
    hlines = nothing,
)
    # Convert nothing vlines/hlines to PrettyTables-friendly defaults
    opts = (
        show_row_number = show_row_number,
        crop = crop,
        show_header = show_header,
        alignment = alignment,
    )
    if vlines === nothing && hlines === nothing
        pretty_table(io, table; opts...)
    elseif vlines === nothing
        pretty_table(io, table; opts..., hlines = hlines)
    elseif hlines === nothing
        pretty_table(io, table; opts..., vlines = vlines)
    else
        pretty_table(io, table; opts..., vlines = vlines, hlines = hlines)
    end
    nothing
end

end # module
