from synbio.polymers import DNA


def seqrecord_to_DNA(record):
    """
    """
    # TODO: docstring
    dna_obj = DNA(record.seq)
    annotation = Part(seq=dna_obj, name=record.name, kind='source')
    return dna_obj


def seqfeature_to_Part(feature):
    """
    """
    # TODO: docstring
    # get location
    start = feature.location.start
    end = feature.location.end
    strand = "FWD" if feature.location.strand == 1 else "REV"
    location = Location(start, end, strand)
    # get modifiers
    name = feature.qualifiers.get('label', '???')
    kind = feature.type
    metadata = feature.qualifiers
    return partial(Part, name=name, kind=kind, location=location, metadata=metadata)


def get_features_by_key(record, key="gene"):
    """
    A function that filters Biopython SeqRecord.features by key
    """
    # TODO: docstring

    def feat_name(feat):
        return feat.qualifiers['locus_tag'][0]

    def feat_seq(feat, record):
        return str(feat.extract(record).seq)

    # return dict of name: sequence for each feature matching key
    return {
        gene_name(feat): gene_seq(feat, record)
        for feat in record.features if feat.type == key
    }
