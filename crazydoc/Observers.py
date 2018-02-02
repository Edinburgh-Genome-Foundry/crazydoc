from .conf import conf
from .biotools import sequence_to_record, sequence_to_annotated_record

class StyleObserver:

    def __init__(self):
        pass

    def process_record_features(self, record):
        for feature in record.features:
            self.process_feature(feature)

    def process_feature(self, feature):
        if self.name not in feature.qualifiers:
            return
        value = feature.qualifiers[self.name]
        label = ''
        if 'label' in feature.qualifiers:
            label = feature.qualifiers['label'] + '; '
        label += self.name
        if value != True:
            label += ": " + str(value)
        feature.qualifiers['label'] = label

    def msword_runs_to_record(self, runs):
        features = [[None, '']]
        for run in runs:
            value = self.evaluate(run)

            if value == features[-1][0]:
                features[-1][1] += run.text
            else:
                features.append([value, run.text])
        feature_records = [
            (
                sequence_to_annotated_record(text, **{self.name: value_})
                if value_
                else sequence_to_record(text)
            )
            for (value_, text) in features
        ]
        record = feature_records[0]
        if len(feature_records) > 1:
            record = sum(feature_records[1:], record)
        return record


class ColorObserver(StyleObserver):

    def process_feature(self, feature):
        if self.name not in feature.qualifiers:
            return
        color = feature.qualifiers[self.name]
        feature.qualifiers['color'] = color
        feature.qualifiers['ApEinfo_revcolor'] = color
        feature.qualifiers['ApEinfo_fwdcolor'] = color

class Italic(StyleObserver):
    name = 'italic'
    def evaluate(self, run):
        return run.italic


class Bold(StyleObserver):
    name = 'bold'
    def evaluate(self, run):
        return run.bold

class Underline(StyleObserver):
    name = 'underline'
    def evaluate(self, run):
        return run.underline


class FontColor(ColorObserver):
    name = 'font_color'
    def evaluate(self, run):
        color = str(run.font.color.rgb)
        if color in ['None', '000000']:
            return False
        else:
            return "#" + color


class HighlightColor(ColorObserver):
    name = 'highlight_color'

    def evaluate(self, run):
        color = run.font.highlight_color
        if color is None:
            return False
        else:
            return conf['color_theme'][color._member_name]


class UpperCase(StyleObserver):
    name = 'upper_case'

    def evaluate(self, run):
        return (run.text == run.text.upper())


class LowerCase(StyleObserver):
    name = 'lower_case'

    def evaluate(self, run):
        return (run.text == run.text.upper())
