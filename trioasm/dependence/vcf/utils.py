"""
Utilities for VCF files.
"""

def walk_together(*readers, **kwargs):
    """
    Simultaneously iteratate over two or more VCF readers. For each 
    genomic position with a variant, return a list of size equal to the number 
    of VCF readers. This list contains the VCF record from readers that have
    this variant, and None for readers that don't have it. 
    The caller must make sure that inputs are sorted in the same way and use the 
    same reference otherwise behaviour is undefined.

    Args:
        vcf_record_sort_key: function that takes a VCF record and returns a 
            tuple that can be used as a key for comparing and sorting VCF 
            records across all readers. This tuple defines what it means for two 
            variants to be equal (eg. whether it's only their position or also 
            their allele values), and implicitly determines the chromosome 
            ordering since the tuple's 1st element is typically the chromosome 
            name (or calculated from it).
    """
    if 'vcf_record_sort_key' in kwargs:
        get_key = kwargs['vcf_record_sort_key']
    else:
        get_key = lambda r: (r.CHROM, r.POS) #, r.REF, r.ALT)

    nexts = []
    for reader in readers:
        try:
            nexts.append(next(reader))
        except StopIteration:
            nexts.append(None)

    min_k = (None,)   # keep track of the previous min key's contig
    while any([r is not None for r in nexts]):
        next_idx_to_k = dict(
            (i, get_key(r)) for i, r in enumerate(nexts) if r is not None)
        keys_with_prev_contig = [
            k for k in list(next_idx_to_k.values()) if k[0] == min_k[0]]

        if any(keys_with_prev_contig):
            min_k = min(keys_with_prev_contig)   # finish previous contig
        else:
            min_k = min(next_idx_to_k.values())   # move on to next contig

        min_k_idxs = set([i for i, k in list(next_idx_to_k.items()) if k == min_k])
        yield [nexts[i] if i in min_k_idxs else None for i in range(len(nexts))]

        for i in min_k_idxs:
            try:
                nexts[i] = next(readers[i])
            except StopIteration:
                nexts[i] = None


def trim_common_suffix(*sequences):
    """
    Trim a list of sequences by removing the longest common suffix while
    leaving all of them at least one character in length.

    Standard convention with VCF is to place an indel at the left-most
    position, but some tools add additional context to the right of the
    sequences (e.g. samtools). These common suffixes are undesirable when
    comparing variants, for example in variant databases.

        >>> trim_common_suffix('TATATATA', 'TATATA')
        ['TAT', 'T']

        >>> trim_common_suffix('ACCCCC', 'ACCCCCCCC', 'ACCCCCCC', 'ACCCCCCCCC')
        ['A', 'ACCC', 'ACC', 'ACCCC']

    """
    if not sequences:
        return []
    reverses = [seq[::-1] for seq in sequences]
    rev_min = min(reverses)
    rev_max = max(reverses)
    if len(rev_min) < 2:
        return sequences
    for i, c in enumerate(rev_min[:-1]):
        if c != rev_max[i]:
            if i == 0:
                return sequences
            return [seq[:-i] for seq in sequences]
    return [seq[:-(i + 1)] for seq in sequences]
