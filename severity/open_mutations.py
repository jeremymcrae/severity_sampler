

MISSENSE_CQ = set(["missense_variant", "stop_lost", "inframe_deletion",
    "inframe_insertion", "coding_sequence_variant", "protein_altering_variant"])

LOF_CQ = set(["stop_gained", "splice_acceptor_variant",
    "splice_donor_variant", "frameshift_variant", "initiator_codon_variant",
    "start_lost", "conserved_exon_terminus_variant"])

def open_mutations(path, indels=False):
    
    genes = {}
    with open(path) as handle:
        header = handle.readline().strip().split('\t')
        chrom_idx = header.index('chrom')
        pos_idx = header.index('pos')
        ref_idx = header.index('ref')
        alt_idx = header.index('alt')
        symbol_idx = header.index('symbol')
        cq_idx = header.index('consequence')
        
        for line in handle:
            line = line.strip().split('\t')
            
            chrom, pos, ref, alt, symbol, cq = line[chrom_idx], \
                line[pos_idx], line[ref_idx], line[alt_idx], \
                line[symbol_idx], line[cq_idx]
            
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
