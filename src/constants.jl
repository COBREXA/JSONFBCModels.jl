
"""
Various magic values and identifiers that are often used in JSON models.
"""
module constants
default_reaction_bound = 1e3

keynames = (
    reactions = ["reactions", "rxns", "RXNS", "REACTIONS", "Reactions", "Rxns"],
    metabolites = ["metabolites", "mets", "METS", "METABOLITES", "Metabolites", "Mets"],
    genes = ["genes", "GENES", "Genes"],
)
end
