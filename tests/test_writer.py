"""
Tests for Crazydoc writer
"""
import os
import pytest
from crazydoc import CrazydocParser, write_crazydoc
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation

example_path = os.path.join("tests", "data", "test_sample.docx")


def test_writer_doesnt_change_feature_locations(tmpdir):
    parser = CrazydocParser(
        ["highlight_color", "bold", "underline", "font_color", "lower_case"]
    )
    seqs = parser.parse_doc_file(example_path)
    for r in seqs:
        # Each feature must have a unique ID
        for i, f in enumerate(r.features):
            f.qualifiers["uID"] = f"feature{i}"
    output_path = os.path.join(
        str(tmpdir), "test_writer_doesnt_change_feature_locations.docx"
    )
    # Write the test sequences
    write_crazydoc(seqs, "uID", output_path)
    # Re-read the test sequences
    written_seqs = parser.parse_doc_file(output_path)
    # Check the locations are unchanged
    for seq, written_seq in zip(seqs, written_seqs):
        locations = [(int(f.location.start), int(f.location.end)) for f in seq.features]
        locations.sort()
        written_locations = [
            (int(f.location.start), int(f.location.end)) for f in written_seq.features
        ]
        written_locations.sort()
        assert locations == written_locations


########################################################################################
# Custom seq example
# too many groups
r = SeqRecord(Seq("abcdefghijklmnop"), id="eg")
for i in range(6):
    f = SeqFeature(SimpleLocation(0, 2), qualifiers={"note": i})
    r.features.append(f)


def test_too_many_groups(tmpdir):
    output_path = os.path.join(str(tmpdir), "test_too_many_groups.docx")
    with pytest.raises(AssertionError):
        write_crazydoc(r, "note", output_path)


########################################################################################
# custom seq example
# too many features in a group
r = SeqRecord(Seq("abcdefghijklmnop"), id="eg")
for i in range(15):
    f = SeqFeature(SimpleLocation(i, i + 1), qualifiers={"note": i})
    r.features.append(f)


def test_too_many_features(tmpdir):
    output_path = os.path.join(str(tmpdir), "test_too_many_features.docx")
    with pytest.raises(AssertionError):
        write_crazydoc(r, "note", output_path)
