from functools import partial

from synbio.polymers import DNA
from synbio.annotations import Part, Location


def seqrecord_to_DNA(record):
    """
    """
    # TODO: docstring
    dna_obj = DNA(record.seq)
    _ = Part(seq=dna_obj, name=record.name, kind='source')
    _ = [
        seqfeature_to_Part(feat)(seq=dna_obj)
        for feat in record.features
    ]

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
    if isinstance(name, list):
        name = name[0]
    kind = feature.type
    metadata = feature.qualifiers
    return partial(Part, name=name, kind=kind,
                   location=location, metadata=metadata)


def get_features_by_key(record, key="gene"):
    """
    A function that filters Biopython SeqRecord.features by key
    """
    # TODO: docstring

    def feat_name(feat):
        return feat.qualifiers['locus_tag'][0]

    def feat_seq(feat, rec):
        return str(feat.extract(rec).seq)

    # return dict of name: sequence for each feature matching key
    return {
        feat_name(feat): feat_seq(feat, record)
        for feat in record.features
        if feat.type == key
    }
