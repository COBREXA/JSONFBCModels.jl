
import PikaParser as PP

"""
`PikaParser.jl` grammar for stringy GRR expressions.
"""
const grr_grammar = begin
    # characters that typically form the identifiers
    isident(x::Char) =
        isletter(x) ||
        isdigit(x) ||
        x == '_' ||
        x == '-' ||
        x == ':' ||
        x == '.' ||
        x == '@' ||
        x == '#' ||
        x == '\'' ||
        x == '[' ||
        x == ']'

    # scanner helpers
    eat(p) = m -> begin
        last = 0
        for i in eachindex(m)
            p(m[i]) || break
            last = i
        end
        return last
    end

    # eat one of keywords
    kws(w...) = m -> begin
        last = eat(isident)(m)
        m[begin:last] in w ? last : 0
    end

    PP.make_grammar(
        [:expr],
        PP.flatten(
            Dict(
                :space => PP.first(PP.scan(eat(isspace)), PP.epsilon),
                :id => PP.scan(eat(isident)),
                :orop =>
                    PP.first(PP.tokens("||"), PP.token('|'), PP.scan(kws("OR", "or"))),
                :andop => PP.first(
                    PP.tokens("&&"),
                    PP.token('&'),
                    PP.scan(kws("AND", "and")),
                ),
                :expr => PP.seq(:space, :orexpr, :space, PP.end_of_input),
                :orexpr => PP.first(
                    :or => PP.seq(:andexpr, :space, :orop, :space, :orexpr),
                    :andexpr,
                ),
                :andexpr => PP.first(
                    :and => PP.seq(:baseexpr, :space, :andop, :space, :andexpr),
                    :baseexpr,
                ),
                :baseexpr => PP.first(
                    :id,
                    :parenexpr => PP.seq(
                        PP.token('('),
                        :space,
                        :orexpr,
                        :space,
                        PP.token(')'),
                    ),
                ),
            ),
            Char,
        ),
    )
end

grr_grammar_open(m, _) =
    m.rule == :expr ? Bool[0, 1, 0, 0] :
    m.rule == :parenexpr ? Bool[0, 0, 1, 0, 0] :
    m.rule in [:or, :and] ? Bool[1, 0, 0, 0, 1] :
    m.rule in [:andexpr, :orexpr, :notexpr, :baseexpr] ? Bool[1] :
    (false for _ in m.submatches)

grr_grammar_fold(m, _, subvals) =
    m.rule == :id ? Expr(:call, :gene, String(m.view)) :
    m.rule == :and ? Expr(:call, :and, subvals[1], subvals[5]) :
    m.rule == :or ? Expr(:call, :or, subvals[1], subvals[5]) :
    m.rule == :parenexpr ? subvals[3] :
    m.rule == :expr ? subvals[2] : isempty(subvals) ? nothing : subvals[1]

"""
$(TYPEDSIGNATURES)

Parses a JSON-ish data reference to a `Expr`-typed gene association. Contains
"calls" to `gene`, `and` and `or` functions that describe the association.
"""
function parse_gene_association(str::String)::Maybe{Expr}
    all(isspace, str) && return nothing
    tree = PP.parse_lex(grr_grammar, str)
    match = PP.find_match_at!(tree, :expr, 1)
    match > 0 || throw(DomainError(str, "cannot parse GRR"))
    PP.traverse_match(tree, match, open = grr_grammar_open, fold = grr_grammar_fold)
end

"""
$(TYPEDSIGNATURES)

Evaluate the gene association expression with the reference values given by the
`val` function.
"""
function eval_gene_association(ga::Expr, val::Function)::Bool
    (ga.head == :call && length(ga.args) >= 2) ||
        throw(DomainError(ga, "invalid gene association expr"))
    if ga.args[1] == :gene && length(ga.args) == 2
        val(ga.args[2])
    elseif ga.args[1] == :and
        all(eval_gene_association.(ga.args[2:end], Ref(val)))
    elseif ga.args[1] == :or
        any(eval_gene_association.(ga.args[2:end], Ref(val)))
    else
        throw(DomainError(ga, "unsupported gene association function"))
    end
end

"""
$(TYPEDSIGNATURES)

A helper for producing predictable unique sequences. Might be faster if
compacting would be done directly in sort().
"""
function sortunique(x)
    o = collect(x)
    sort!(o)
    put = prevind(o, firstindex(o))
    for i in eachindex(o)
        if put >= firstindex(o) && o[i] == o[put]
            # we already have this one
            continue
        else
            put = nextind(o, put)
            if put != i
                o[put] = o[i]
            end
        end
    end
    o[begin:put]
end

"""
$(TYPEDSIGNATURES)

Convert the given gene association expression to DNF.
"""
function flatten_gene_association(ga::Expr)::A.GeneAssociationDNF
    function fold_and(dnfs::Vector{Vector{Vector{String}}})::Vector{Vector{String}}
        if isempty(dnfs)
            [String[]]
        else
            sortunique(
                sortunique(String[l; r]) for l in dnfs[1] for r in fold_and(dnfs[2:end])
            )
        end
    end

    (ga.head == :call && length(ga.args) >= 2) ||
        throw(DomainError(ga, "invalid gene association expr"))
    if ga.args[1] == :gene && length(ga.args) == 2
        [[ga.args[2]]]
    elseif ga.args[1] == :and
        fold_and(flatten_gene_association.(ga.args[2:end]))
    elseif ga.args[1] == :or
        sortunique(vcat(flatten_gene_association.(ga.args[2:end])...))
    else
        throw(DomainError(ga, "unsupported gene association function"))
    end
end

"""
$(TYPEDSIGNATURES)

Formats a DNF gene association as a `String`.
"""
function format_gene_association_dnf(
    grr::A.GeneAssociationDNF;
    and = " && ",
    or = " || ",
)::String
    return join(("(" * join(gr, and) * ")" for gr in grr), or)
end
