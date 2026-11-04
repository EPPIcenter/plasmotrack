"""Values of an Analysis taken over its Chains: its Infections, and per-Chain means pooled over every Draw."""


def infection_ids(analysis):
    """The Infection ids of every Chain of `analysis`, as a set."""
    return {infection for store in analysis.runs.values()
            for (infection,) in store.structure("infections").project("infection").fetchall()}


def pooled_mean(chains, key, draws):
    """The mean over every Chain's Draws of a per-Chain mean: `chains` is [{chain, key: mean}], `draws` maps each
    Chain to its number of Draws, and each Chain weighs by its Draws. None without Draws."""
    total = sum(draws[c["chain"]] for c in chains)
    return sum(c[key] * draws[c["chain"]] for c in chains) / total if total else None
