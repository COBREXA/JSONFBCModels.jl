
"""
A named tuple that contains the magic values that are used globally for
whatever purposes.
"""
const constants = (
    default_reaction_bound = 1e3,
    keynames = (
        reactions = ["reactions", "rxns", "RXNS", "REACTIONS", "Reactions", "Rxns"],
        metabolites = ["metabolites", "mets", "METS", "METABOLITES", "Metabolites", "Mets"],
        genes = ["genes", "GENES", "Genes"],
    ),
)

const Maybe = A.Maybe
