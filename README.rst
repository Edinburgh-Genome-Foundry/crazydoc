.. raw:: html

    <p align="center">
    <img alt="crazydoc Logo" title="crazydoc Logo" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/crazydoc/master/docs/title.png" width="550">
    <br /><br />
    </p>

.. image:: https://travis-ci.org/Edinburgh-Genome-Foundry/crazydoc.svg?branch=master
   :target: https://travis-ci.org/Edinburgh-Genome-Foundry/crazydoc
   :alt: Travis CI build status

.. image:: https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/crazydoc/badge.svg?branch=master
   :target: https://coveralls.io/github/Edinburgh-Genome-Foundry/crazydoc?branch=master


Crazydoc is a Python library to parse one of the most common DNA representation formats: the joyfully coloured and stylishly annotated MS-Word document.

.. raw:: html

    <p align="center">
    <img src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/crazydoc/master/docs/screenshot.png" width="600">
    </p>

While other standards such as FASTA or Genbank are better supported by modern sequence editors, none enjoys the same popularity among molecular biologist as MS-Word's ``.docx`` format, which is limited only by the sophistication and creativity of the user.

Relying on a loose syntax and unclear specifications, this format has however suffered from a lack of support in the developers community and is generally incompatible with mainstream software pipelines. This library allows to convert MS-Word DNA sequences to more computing friendly formats: Biopython records, FASTA, or annotated Genbanks.

Usage
-----

To obtain all sequences contained in a docx as annotated Biopython records (such as `this one <https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/crazydoc/master/examples/example.docx>`_):

.. code:: python

    from crazydoc import CrazydocParser
    parser = CrazydocParser(['highlight_color', 'bold', 'underline'])
    biopython_records = parser.parse_doc_file("./example.docx")

You can then plot the obtained records:

.. code:: python

    from crazydoc import CrazydocSketcher
    sketcher = CrazydocSketcher()
    for record in biopython_records:
        sketch = sketcher.translate_record(record)
        ax, _ = sketch.plot()
        ax.set_title(record.id)
        ax.figure.savefig('%s.png' % record.id)

.. raw:: html

    <p align="center">
    <img src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/crazydoc/master/docs/records_plots.png" width="800">
    </p>

To write the sequences down as Genbank records, with annotations:

.. code:: python

    from crazydoc import records_to_genbank
    records_to_genbank(biopython_records)

Installation
-------------

(soon) You can install crazydoc through PIP

.. code::

    sudo pip install crazydoc

Alternatively, you can unzip the sources in a folder and type

.. code::

    sudo python setup.py install

License = MIT
--------------

Crazydoc is an open-source software originally written at the `Edinburgh Genome Foundry <http://genomefoundry.org>`_ by `Zulko <https://github.com/Zulko>`_ and `released on Github <https://github.com/Edinburgh-Genome-Foundry/crazydoc>`_ under the MIT licence (¢ Edinburg Genome Foundry).

Everyone is welcome to contribute !

More biology software
---------------------

.. image:: https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Edinburgh-Genome-Foundry.github.io/master/static/imgs/logos/egf-codon-horizontal.png
  :target: https://edinburgh-genome-foundry.github.io/

Crazydoc is part of the `EGF Codons <https://edinburgh-genome-foundry.github.io/>`_ synthetic biology software suite for DNA design, manufacturing and validation.
