def get_features_by_key(record, key="gene"):
    '''
    A function that filters Biopython SeqRecord.features by key
    '''
    def feat_name(feat):
        return feat.qualifiers['locus_tag'][0]
    def feat_seq(feat, record):
        return str(feat.extract(record).seq)
    # return dict of name: sequence for each feature matching key
    return {
        gene_name(feat): gene_seq(feat, record)
            for feat in record.features if feat.type == key
    }
