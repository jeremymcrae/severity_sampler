

MISSENSE_CQ = set(["missense_variant", "stop_lost", "inframe_deletion",
    "inframe_insertion", "coding_sequence_variant", "protein_altering_variant"])

LOF_CQ = set(["stop_gained", "splice_acceptor_variant",
    "splice_donor_variant", "frameshift_variant", "initiator_codon_variant",
    "start_lost", "conserved_exon_terminus_variant"])

def open_mutations(path, indels=False):
    
    genes = {}
    with open(path) as handle:
        header = handle.readline().strip().split('\t')
        
        for line in handle:
            _, chrom, pos, ref, alt, symbol, cq, *_ = line.strip().split('\t')
            
            # ignore noncoding
            if cq not in MISSENSE_CQ | LOF_CQ:
                continue
            
            if symbol not in genes:
                genes[symbol] = []
            
            if not indels and (len(ref) > 1 or len(alt) > 1):
                continue
            
            genes[symbol].append({'chrom': chrom, 'position': int(pos), 'ref': ref,
                'alt': alt, 'consequence': cq})
    
    return genes
