import os

from copy import deepcopy
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

try:
    # Biopython <1.78
    from Bio.Alphabet import DNAAlphabet
    from Bio.Alphabet import generic_protein

    has_dna_alphabet = True
except ImportError:
    # Biopython >=1.78
    has_dna_alphabet = False
from Bio.SeqFeature import SeqFeature, FeatureLocation


def annotate_record(
    seqrecord, location="full", feature_type="misc_feature", margin=0, **qualifiers
):
    """Add a feature to a Biopython SeqRecord.

    Parameters
    ----------

    seqrecord
      The Biopython SeqRecord to be annotated.

    location
      Either (start, end) or (start, end, strand). (strand defaults to +1).

    feature_type
      The type associated with the feature.

    margin
      Number of extra bases added on each side of the given location.

    qualifiers
      Dictionary that will be the Biopython feature's `qualifiers` attribute.
    """
    if location == "full":
        location = (margin, len(seqrecord) - margin)

    strand = location[2] if len(location) == 3 else 1
    seqrecord.features.append(
        SeqFeature(
            FeatureLocation(location[0], location[1], strand),
            qualifiers=qualifiers,
            type=feature_type,
        )
    )


def sequence_to_record(
    sequence, id="<unknown id>", name="<unknown name>", is_protein=False
):
    """Return a SeqRecord of the sequence, ready to be Genbanked."""
    if has_dna_alphabet:  # Biopython <1.78
        if is_protein:
            seq = Seq(sequence, alphabet=generic_protein)
        else:
            seq = Seq(sequence, alphabet=DNAAlphabet())
    else:
        seq = Seq(sequence)

    seqrecord = SeqRecord(seq, id=id, name=name)
    if is_protein:
        seqrecord.annotations["molecule_type"] = "protein"
    else:
        seqrecord.annotations["molecule_type"] = "DNA"

    return seqrecord


def sequence_to_annotated_record(
    sequence,
    id="<unknown id>",
    name="<unknown name>",
    feature_type="misc_feature",
    is_protein=False,
    **qualifiers
):
    record = sequence_to_record(sequence, is_protein=is_protein)
    annotate_record(record, feature_type=feature_type, **qualifiers)
    return record


genetic_code = "ABCDGHKMNRSTVWY"
nucleotides_set = set(genetic_code + genetic_code.lower())
protein_code = "ACDEFGHIKLMNPQRSTVWYX"
aminoacids_set = set(protein_code + protein_code.lower())


def string_is_sequence(string, is_protein=False):
    if is_protein:
        return len(string) and (set(string) <= aminoacids_set)
    else:
        return len(string) and (set(string) <= nucleotides_set)


def records_to_genbank(records, path=".", extension=None, is_protein=False):
    if extension is None:
        if is_protein:
            extension = ".gp"
        else:
            extension = ".gbk"
    name_duplicates = {}
    for record in records:
        name = record.id.replace("/", "-")
        if name in name_duplicates:
            name_duplicates[name] += 1
            name = "%s_%03d" % (name, name_duplicates[name])
        out_record = deepcopy(record)
        out_record.name = record.name[:20]
        SeqIO.write(out_record, os.path.join(path, name + extension), "genbank")
