try:
    from dna_features_viewer import BiopythonTranslator

    class CrazydocSketcher(BiopythonTranslator):
        """Custom sequence plotter for records from crazydoc.

        Examples:
        ---------
        >>> # plot and export the schema of records from crazydoc
        >>> translator = CrazydocGraphicTranslator()
        >>> for record in records:
        >>>     ax, _ = translator.translate_record(record).plot();
        >>>     ax.figure.savefig(record.id + ".png")
        """

        default_feature_color = "#ffffff"

        def compute_feature_label(self, feature):
            return feature.qualifiers.get("label", None)


except ImportError:

    class CrazydocSketcher:
        """Unavailable class: install DnaFeaturesViewer.

        Try: pip install dna_features_viewer
        """

        def __init__(self):
            raise ImportError(
                "The graphic translator requires "
                "DnaFeaturesViewer to be installed.\n"
                "(try: pip install dna_features_viewer)."
            )
