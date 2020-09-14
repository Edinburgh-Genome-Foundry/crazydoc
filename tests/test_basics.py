import os
import matplotlib

matplotlib.use("Agg")
from crazydoc import CrazydocParser, CrazydocSketcher, records_to_genbank


example_path = os.path.join("tests", "data", "test_sample.docx")


def test_crazydocparser(tmpdir):
    parser = CrazydocParser(
        ["highlight_color", "bold", "underline", "lower_case", "font_color"]
    )
    biopython_records = parser.parse_doc_file(example_path)
    record_features = [
        [(f.location.start, f.location.end, f.qualifiers) for f in record.features]
        for record in biopython_records
    ]
    assert len(
        [
            start
            for start, end, qual in record_features[0]
            if (start, end, qual.get("lower_case", False)) == (1237, 1253, True)
        ]
    )
    records_to_genbank(biopython_records, path=str(tmpdir))

    sketcher = CrazydocSketcher()
    for record in biopython_records:
        sketch = sketcher.translate_record(record)
        ax, _ = sketch.plot()
        ax.set_title(record.id)
        filepath = os.path.join(str(tmpdir), "%s.png" % record.id)
        ax.figure.savefig(filepath, bbox_inches="tight")
