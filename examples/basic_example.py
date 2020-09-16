from crazydoc import CrazydocParser, CrazydocSketcher, records_to_genbank
import os

parser = CrazydocParser(["highlight_color", "bold", "underline"])
biopython_records = parser.parse_doc_file("example.docx")

output_dir = "examples_outputs"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
records_to_genbank(biopython_records, path=output_dir)

sketcher = CrazydocSketcher()
for record in biopython_records:
    sketch = sketcher.translate_record(record)
    ax, _ = sketch.plot()
    ax.set_title(record.id)
    filepath = os.path.join("examples_outputs", "%s.png" % record.id)
    ax.figure.savefig(filepath, bbox_inches="tight")

print("The genbanks and images files are now in folder examples_ouputs/")
