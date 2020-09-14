from .conf import conf
from .biotools import sequence_to_record, sequence_to_annotated_record


class StyleObserver:
    """Generic class for observing style-based annotations in sequences.

    The provided subclasses each observe one particular type of DNA sequence
    annotation, such as the highlight color, bold text, underlines, etc.
    """

    def __init__(self):
        pass

    def process_record_features(self, record):
        for feature in record.features:
            self.process_feature(feature)

    def process_feature(self, feature):
        if self.name not in feature.qualifiers:
            return
        value = feature.qualifiers[self.name]
        label = ""
        if "label" in feature.qualifiers:
            label = feature.qualifiers["label"] + "; "
        label += self.name
        if not isinstance(value, bool):
            label += ": " + str(value)
        feature.qualifiers["label"] = label

    def aggregate_features_from_runs(self, runs):
        features = [[None, ""]]
        for run in runs:
            value = self.evaluate(run)
            text = run.text.replace(" ", "")
            if value == features[-1][0]:
                features[-1][1] += text
            else:
                features.append([value, text])
        return features

    def msword_runs_to_record(self, runs, is_protein=False):
        feature_records = [
            (
                sequence_to_annotated_record(
                    text, is_protein=is_protein, **{self.name: val}
                )
                if val
                else sequence_to_record(text, is_protein=is_protein)
            )
            for (val, text) in self.aggregate_features_from_runs(runs)
        ]
        record = feature_records[0]
        if len(feature_records) > 1:
            record = sum(feature_records[1:], record)
        return record


class ColorObserver(StyleObserver):
    """Subclass for color observers."""

    def process_feature(self, feature):
        if self.name not in feature.qualifiers:
            return
        color = feature.qualifiers[self.name]
        for field in ["color", "ApEinfo_revcolor", "ApEinfo_fwdcolor"]:
            feature.qualifiers[field] = color


class CharactersObserver(StyleObserver):
    """Subclass for character-by-character observers."""

    def aggregate_features_from_runs(self, runs):
        features = [[None, ""]]
        text = "".join([r.text for r in runs])
        for character in text:
            value = self.evaluate(character)
            if value == features[-1][0]:
                features[-1][1] += character
            else:
                features.append([value, character])
        return features


class Italic(StyleObserver):
    """Captures italic text."""

    name = "italic"

    def evaluate(self, run):
        """Return whether the run has italic style"""
        return run.italic


class Bold(StyleObserver):
    """Captures bold text."""

    name = "bold"

    def evaluate(self, run):
        """Return whether the run has bold style"""
        return run.bold


class Underline(StyleObserver):
    """Captures underlined text."""

    name = "underline"

    def evaluate(self, run):
        """Return whether the run has underline style"""
        return run.underline


class FontColor(ColorObserver):
    """Captures text with non-black font color."""

    name = "font_color"

    def evaluate(self, run):
        """Return False if no color, else the #ae60bf color."""
        color = str(run.font.color.rgb)
        if color in ["None", "000000"]:
            return False
        else:
            return "#" + color


class HighlightColor(ColorObserver):
    """Captures text with a background-highlighting color."""

    name = "highlight_color"

    def evaluate(self, run):
        """Return False if no background color, else the #ae60bf color."""
        color = run.font.highlight_color
        if color is None:
            return False
        else:
            return conf["color_theme"][color._member_name]


class UpperCase(CharactersObserver):
    """Captures upper-case text."""

    name = "upper_case"

    def evaluate(self, character):
        """Return whether the character is upper"""
        return character == character.upper()


class LowerCase(CharactersObserver):
    """Captures lower-case text."""

    name = "lower_case"

    def evaluate(self, character):  #
        """Return whether the character is lower"""
        return character == character.lower()
