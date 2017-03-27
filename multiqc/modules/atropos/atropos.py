"""MultiQC module that reads JSON output from Atropos.

Note: this code is mostly borrowed from the Cutadapt [1] and FastQC [2] MultiQC
modules.
1. https://github.com/ewels/MultiQC/blob/master/multiqc/modules/fastqc/fastqc.py
2. https://github.com/ewels/MultiQC/blob/master/multiqc/modules/cutadapt/cutadapt.py
"""
from __future__ import print_function, division, absolute_import
from collections import OrderedDict
import logging
import io
import math
import operator
import os
import json
import re
import scipy.stats

from multiqc import config
from multiqc.plots import linegraph, bargraph
from multiqc.modules.base_module import BaseMultiqcModule

## Hard-coded URLs
# TODO: eventually these need pointers to specific sections of the docs

ATROPOS_GITHUB_URL = "https://github.com/jdidion/atropos"
ATROPOS_DOC_URL = "http://atropos.readthedocs.org/en/latest/guide.html"

## Assets

ATROPOS_CSS = {
    'assets/css/multiqc_atropos.css' :
    os.path.join(os.path.dirname(__file__), 'assets', 'css', 'multiqc_fastqc.css')
}
ATROPOS_JS = {
    'assets/js/multiqc_atropos.js' :
    os.path.join(os.path.dirname(__file__), 'assets', 'js', 'multiqc_fastqc.js')
}
ATROPOS_COLORS = {
    'pass': '#5cb85c', 'warn': '#f0ad4e', 'fail': '#d9534f', 'default': '#999'
}

ATROPOS_TRIMMED_LENGTH = """
<p>This plot shows the number of reads with certain lengths of adapter trimmed.
Obs/Exp shows the raw counts divided by the number expected due to sequencing
errors. A defined peak may be related to adapter length. See the
<a href="{}" target="_blank">Atropos documentation</a> for more information on
how these numbers are generated.</p>""".format(ATROPOS_DOC_URL)

ATROPOS_PASSFAILS = "<script type="text/javascript">atropos_passfails = {};</script>"

# Initialise the logger
log = logging.getLogger(__name__)

## Module

class MultiqcModule(BaseMultiqcModule):
    """Atropos module class. Loads JSON ouptut for read trimming stats (which
    are equivalent to Cutadapt) and pre- and post-trimming stats (which are
    very similar to FastQC).
    
    Args:
        from_config: Whether to load data from configured locations. This should
            only be set to False in tests.
    """
    def __init__(self, from_config=True):
        # Initialise the parent object
        super().__init__(
            name='Atropos', anchor='atropos',
            href=ATROPOS_GITHUB_URL,
            info="is a general-purpose NGS pre-processing tool that"\
                 "specializes in adatper- and quality-trimming.")
        
        # List of sample IDs
        self.atropos_sample_ids = []
        # Collections of data dicts
        self.atropos_data = OrderedDict()
        self.atropos_trim_data = OrderedDict()
        self.atropos_qc_data = dict(pre=OrderedDict(), post=OrderedDict())
        # Section objects
        self.sections = [
            PerBaseQuality(),
            #PerTileQuality()
            PerSequenceQuality(),
            #PerBaseContent(),
            PerSequenceGC(),
            PerBaseN(),
            #SequenceLength(),
            #Duplication(),
            #Contaminants(),
            #Kmers()
        ]
        # Sections of the report
        self.section_html = []
        # Add to self.css and self.js to be included in template
        self.css = ATROPOS_CSS
        self.js = ATROPOS_JS
        # Colours to be used for plotting lines
        self.status_colours = ATROPOS_JS
            
        if from_config:
            self.init_from_config()
            self.atropos_report()
    
    def init_from_config(self):
        """Initialize module from JSON files discovered via MultiQC
        configuration.
        """
        for file_dict in self.find_log_files(
                config.sp['atropos'], patterns=dict(fn='*.json'),
                filehandles=True):
            fileobj = file_dict['f']
            data = json.load(fileobj)
            self.add_atropos_data(data, fileobj)
        
        num_samples = len(self.atropos_sample_ids)
        if num_samples == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning
        else:
            log.info("Found {} reports".format(num_samples))
    
    def add_atropos_data(self, data=None, data_source=None):
        """For each sample, Atropos generates a <sample>.json file. That
        file has summary information (numbers of files processed, description of
        the inputs (e.g. FASTA vs FASTQ, single vs paired-end), and a summary
        of the analyses performed.
        
        Args:
            data: Atropos summary dict.
            data_source: The file from which 'data' was loaded.
        """
        close = False
        if data_source and isinstance(data_source, str):
            data_source = open(data_source)
            close = True
        if data is None:
            if data_source:
                data = json.load(data_source)
                if close:
                    data_source.close()
            else:
                raise ValueError(
                    "One of 'data' or 'data_source' must be provided")
            
        sample_id = data['sample_id']
        self.atropos_sample_ids.append(sample_id)
        self.atropos_data[sample_id] = data
        if data_source:
            self.add_data_source(data_source, sample_id)
        
        if 'trim' in data:
            self.atropos_trim_data[sample_id] = data['trim']
        
        if 'pre' in data:
            for source, pre_data in data['pre'].items():
                # TODO
        
        # TODO: handle multiplexed output
        if 'post' in data and 'NoFilter' in data['post']:
            for source, post_data in data['post']['NoFilter'].items():
                # TODO
        
        # Load qc data
        if 'qc_stats' in data:
            for qc_file in data['qc_stats']:
                with open(qc_file, 'rb') as infile:
                    qc = json.load(infile)
                    qc_id = '{}_{}'.format(sample_id, qc['name'])
                    for phase in ('pre', 'post'):
                        if phase in qc:
                            self.atropos_qc_data[phase][qc_id] = qc[phase]
    
    def atropos_report(self):
        if not self.atropos_sample_ids:
            log.debug("No reports to process; atropos_report raising UserWarning")
            raise UserWarning
        
        # Add to general stats
        if self.atropos_summary_data:
            self.atropos_general()
        
        # Add pre-trim stats
        if self.atropos_qc_data['pre']:
            self.atropos_qc('pre')
        
        # Add trimming stats
        if self.atropos_trim_data:
            self.atropos_trim()
        
        # Add post-trim stats
        if self.atropos_qc_data['post']:
            self.atropos_qc('post')
    
    def atropos_general(self):
        """Add some single-number stats to the basic statistics table at the
        top of the report.
        """
        headers = OrderedDict()
        headers['total_sequences'] = {
            'title': 'M Seqs',
            'description': 'Total Sequences (millions)',
            'min': 0,
            'scale': 'Blues',
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        headers['total_bp'] = {
            'title': 'M bp',
            'description': 'Total Base Pairs (millions)',
            'min': 0,
            'scale': 'Blues',
            'modify': lambda x: x / 1000000,
            'shared_key': 'base_count'
        }
        self.general_stats_addcols(self.atropos_summary_data, headers)
    
    def atropos_qc(self, phase):
        # Add statuses to intro. Note that this is slightly different than
        # FastQC: Atropos reports the relevant statistic, and the threshold
        # for pass/warn/fail is configured in MutliQC (defaulting to
        # the thresholds defined in FastQC).
        statuses = {}
        for section in self.sections:
            section_name = "{}_{}".format(phase, section.name)
            section_data = self.atropos_qc_data[phase]
            statuses[section_name] = section.get_statuses(section_data)
            self.section_html.append(section.plot(section_data))
        
        self.intro += ATROPOS_PASSFAILS.format(json.dumps(statuses))
    
    def atropos_trim(self):
        """ Take the parsed stats from the Cutadapt report and add it to the
        basic stats table at the top of the report.
        """
        headers = {}
        headers['percent_trimmed'] = {
            'title': '% Trimmed',
            'description': '% Total Base Pairs trimmed',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlBu-rev',
            'format': '{:.1f}%'
        }
        header['bp_processed']
        self.general_stats_addcols(self.atropos_data, headers)

        """ Generate the trimming length plot.
        """

        html = ATROPOS_TRIMMED_LENGTH

        pconfig = {
            'id': 'cutadapt_plot',
            'title': 'Lengths of Trimmed Sequences',
            'ylab': 'Counts',
            'xlab': 'Length Trimmed (bp)',
            'xDecimals': False,
            'ymin': 0,
            'tt_label': '<b>{point.x} bp trimmed</b>: {point.y:.0f}',
            'data_labels': [{'name': 'Counts', 'ylab': 'Count'},
                            {'name': 'Obs/Exp', 'ylab': 'Observed / Expected'}]
        }

        html += linegraph.plot(
            [self.cutadapt_length_counts, self.cutadapt_length_obsexp],
            pconfig)

## Utils

class NormalDistribution(object):
    def __init__(self, mean, stdev):
        self.mean = mean
        self.sd2 = sd2 = 2 * self.stdev * self.stdev
    
    def __call__(self, value):
        return (
            math.pow(meth.e, -(math.pow(value - mean, 2) / self.sd2)) /
            math.sqrt(self.sd2 * math.pi))

def ordered_dict(iterable):
    d = OrderedDict()
    for key, value in iterable:
        d[key] = value
    d

## QC sections

class Section(object):
    def get_statuses(self, data, phase):
        return ordered_dict(
            (sample_id, self.get_status(sample_data, phase))
            for sample_id, sample_data in data)
    
    def get_status(self, sample_data, phase):
        if self.threshold_statistic in sample_data:
            stat = sample_data[self.threshold_statistic]
        else:
            stat = self.compute_statistic(sample_data)
        return self.get_status_for(stat, phase)
    
    def get_status_for(self, stat, phase):
        # TODO: Add configuration for thresholds.
        for status, threshold in zip(('fail', 'warn'), self.default_threshold):
            if self.compare(stat, threshold):
                return status
        else:
            return 'pass'
    
    def get_html(self, data):
        plot_data = self.get_plot_data(data)
        return {
            'name': self.display,
            'anchor': self.anchor,
            'content': self.html.format(*self.get_html_variables(data))
        }
    
    def get_html_variables(self, data):
        return self.plot_type.plot(get_plot_data(data), self.plot_config)

## Static HTML/templates

class PerBaseQuality(Section):
    name = 'base_quality'
    display = 'Sequence Quality Histograms'
    anchor = 'atropos_per_base_sequence_quality'
    html = """
<p>The mean quality value across each base position in the read. See the
<a href="{}" target="_bkank">Atropos help</a>.</p>
{{}}""".format(ATROPOS_DOC_URL)
    # TODO: Add boxplots as in FastQC output.
    plot_type = linegraph
    plot_config = {
        'id': 'fastqc_per_base_sequence_quality_plot',
        'title': 'Mean Quality Scores',
        'ylab': 'Phred Score',
        'xlab': 'Position (bp)',
        'ymin': 0,
        'xDecimals': False,
        'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}',
        'colors': self.get_status_cols('per_base_sequence_quality'),
        'yPlotBands': [
            {'from': 28, 'to': 100, 'color': '#c3e6c3'},
            {'from': 20, 'to': 28, 'color': '#e6dcc3'},
            {'from': 0, 'to': 20, 'color': '#e6c3c3'},
        ]
    }

    def get_status_for(self, stat, phase):
        lower_quartile, lowest_median = stat
        if lower_quartile < 5 or lowest_median < 20:
            return 'fail'
        if lower_quartile < 10 or lowest_median < 25:
            return 'warn'
        return 'pass'
    
    def get_plot_data(self, data):
        return ordered_dict(
            (sample_id, hist_to_means(values['base_quality']))
            for sample_id, values in data.items())
    
class PerTileQuality(Section):
    name = 'tile_sequence_quality'
    compare = operator.lt
    threshold_statistic = 'min_tile_versus_mean'
    default_thresholds = (-5, -2)
    
    def plot(self, data):
        # TODO
        pass

class PerSequenceQuality(Section):
    name = 'tile_sequence_quality'
    display = 'Per Sequence Quality Scores'
    anchor = 'atropos_per_sequence_quality_scores'
    compare = operator.lt
    threshold_statistic = 'mode_sequence_quality'
    default_thresholds = (20, 27)
    html = """
<p>The number of reads with average quality scores. Shows if a subset of reads
has poor quality. See the <a href="{}" target="_bkank">Atropos help</a>.</p>
{{}}""".format(ATROPOS_DOC_URL)
    plot_type = linegraph
    plot_config = {
        'id': 'fastqc_per_sequence_quality_scores_plot',
        'title': 'Per Sequence Quality Scores',
        'ylab': 'Count',
        'xlab': 'Mean Sequence Quality (Phred Score)',
        'ymin': 0,
        'xmin': 0,
        'xDecimals': False,
        'colors': self.get_status_cols('per_sequence_quality_scores'),
        'tt_label': '<b>Phred {point.x}</b>: {point.y} reads',
        'xPlotBands': [
            {'from': 28, 'to': 100, 'color': '#c3e6c3'},
            {'from': 20, 'to': 28, 'color': '#e6dcc3'},
            {'from': 0, 'to': 20, 'color': '#e6c3c3'},
        ]
    }
    
    def get_plot_data(self, data):
        return ordered_dict(
            (sample_id, hist_to_means(values['sequence_quality']))
            for sample_id, values in data.items())

class PerBaseContent(Section):
    name = 'bases'
    display = 'Per Base Sequence Content'
    anchor = 'atropos_per_base_sequence_content'
    html = """
<p>The proportion of each base position for which each of the four normal DNA
bases has been called. See the <a href="{}" target="_bkank">FastQC help</a>.</p>
<p class="text-primary"><span class="glyphicon glyphicon-info-sign"></span>
Click a heatmap row to see a line plot for that dataset.</p>
<div id="fastqc_per_base_sequence_content_plot">
    <h5><span class="s_name"><em class="text-muted">rollover for sample name</em></span></h5>
    <div class="fastqc_seq_heatmap_key">
        Position: <span id="fastqc_seq_heatmap_key_pos">-</span>
        <div><span id="fastqc_seq_heatmap_key_t"> %T: <span>-</span></span></div>
        <div><span id="fastqc_seq_heatmap_key_c"> %C: <span>-</span></span></div>
        <div><span id="fastqc_seq_heatmap_key_a"> %A: <span>-</span></span></div>
        <div><span id="fastqc_seq_heatmap_key_g"> %G: <span>-</span></span></div>
    </div>
    <div id="fastqc_seq_heatmap_div" class="fastqc-overlay-plot">
        <div id="fastqc_seq" class="hc-plot">
            <canvas id="fastqc_seq_heatmap" height="100%" width="800px" style="width:100%;"></canvas>
        </div>
    </div>
    <div class="clearfix"></div>
</div>
<script type="text/javascript">
    fastqc_seq_content_data = {{}};
    $(function () {{{{ atropos_seq_content_heatmap(); }}}});
</script>"""
    
    def get_html_variables(self, data):
        """Create the epic HTML for the FastQC sequence content heatmap.
        """
        return json.dumps(ordered_dict(
            (sample_id, values['sequence_gc']['hist'])
            for sample_id, values in data.items()))

class PerSequenceGC(Section):
    name = 'sequence_gc'
    display = 'Per Sequence GC Content'
    anchor = 'atropos_per_sequence_gc_content'
    compare = operator.gt
    threshold_statistic = 'fraction_nonnormal'
    default_thresholds = (0.3, 0.15)
    html =  """
<p>The average GC content of reads. Normal random library typically have a
roughly normal distribution of GC content. See the <a href="{}" target="_bkank">
Atropos help</a>.</p>
<p>The dashed black line shows theoretical GC content: {{}}.</p>
{{}}""".format(ATROPOS_DOC_URL)
    
    def get_html_variables(self, data):
        """Create the HTML for the GC content plot.
        """
        plot_data = ordered_dict(
            (sample_id, values['sequence_gc']['hist'])
            for sample_id, values in data.items())
        
        totals = ordered_dict(
            (sample_id, sum(gc.values()))
            for sample_id, gc in plot_data.items())
        
        def normalize_gc(gc, total):
            return ordered_dict(
                (pct, count * 100 / total)
                for pct, count in gc)
        
        plot_data_norm = ordered_dict(
            (sample_id, normalize_gc(gc, totals[sample_id]))
            for sample_id, gc in plot_data.items())
        
        def theoretical_gc(dist, total):
            nd = NormalDistribution(dist['mean'], dist['stdev'])
            return [nd(i) * total for i in range(101)]
        
        theoretical_data = ordered_dict(
            (sample_id, theoretical_gc(
                values['sequence_gc']['dist'], totals[sequence_id]))
            for sample_id, values in data.items())
        
        plot_config = {
            'id': 'atropos_per_sequence_gc_content_plot',
            'title': 'Per Sequence GC Content',
            'ylab': 'Count',
            'xlab': '%GC',
            'ymin': 0,
            'xmax': 100,
            'xmin': 0,
            'yDecimals': False,
            'tt_label': '<b>{point.x}% GC</b>: {point.y}',
            'colors': self.get_status_cols('per_sequence_gc_content'),
            'data_labels': [
                {'name': 'Percentages', 'ylab': 'Percentage'},
                {'name': 'Counts', 'ylab': 'Count'}
            ]
        }
        
        esconfig = {
            'name': 'Theoretical GC Content',
            'dashStyle': 'Dash',
            'lineWidth': 2,
            'color': '#000000',
            'marker': { 'enabled': False },
            'enableMouseTracking': False,
            'showInLegend': False,
        }
        pconfig['extra_series'] = [ [dict(esconfig)], [dict(esconfig)] ]
        pconfig['extra_series'][0][0]['data'] = theoretical_data
        theoretical_gc_desc = .format(theoretical_gc_name)
        
        return (
            theoretical_gc_name,
            linegraph.plot([plot_data_norm, plot_data], plot_config))

class PerBaseN(Section):
    name = 'per_base_N_content'
    display = 'Per Base N Content'
    anchor = 'fastqc_per_base_n_content'
    compare = operator.gt
    threshold_statistic = 'frac_N'
    default_thresholds = (0.2, 0.05)
    html = """
<p>The percentage of base calls at each position for which an N was called.
See the <a href="{}" target="_bkank">Atropos help</a>.</p>
{{}}""".format(ATROPOS_DOC_URL)
    plot_type = linegraph
    plot_config = {
        'id': 'atropos_per_base_n_content_plot',
        'title': 'Per Base N Content',
        'ylab': 'Percentage N-Count',
        'xlab': 'Position in Read (bp)',
        'yCeiling': 100,
        'yMinRange': 5,
        'ymin': 0,
        'xmin': 0,
        'xDecimals': False,
        'colors': self.get_status_cols('per_base_n_content'),
        'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}%',
        'yPlotBands': [
            {'from': 20, 'to': 100, 'color': '#e6c3c3'},
            {'from': 5, 'to': 20, 'color': '#e6dcc3'},
            {'from': 0, 'to': 5, 'color': '#c3e6c3'},
        ]
    }
    
    def get_plot_data(self, data):
        """Create the HTML for the per base N content plot.
        """
        return ordered_dict(
            (sample_id, values['frac_N'])
            for sample_id, values in data.items())

class SequenceLength(Section):
    name = 'sequence_length'
    display = 'Sequence Length Distribution'
    anchor = 'atropos_sequence_length_distribution'
    threshold_statistic = 'length_range'
    html = """
<p>The distribution of fragment sizes (read lengths) found.
See the <a href="{}" target="_bkank">FastQC help</a>.</p>
{{}}""".format(ATROPOS_DOC_URL)
    all_same_html = """
<p>All samples have sequences of a single length ({} bp).</p>"""
    all_same_within_samples_html = """
<p>All samples have sequences of a single length ({} bp).
See the <a href="#general_stats">General Statistics Table</a>.</p>"""
    
    def get_status_for(self, stat, phase):
        """This always returns 'pass' for post-trimming.
        """
        if phase == 'pre':
            min_len, max_len = stat
            if min_len == 0:
                return 'fail'
            if min_len != max_len:
                return 'warn'
        return 'pass'
    
    def get_html(self, data):
        """Create the HTML for the Sequence Length Distribution plot.
        """
        plot_data = ordered_dict(
            (sample_id, values['sequence_length'])
            for sample_id, values in data.items())
        
        unique_lengths = set()
        multiple_lengths = False
        for lengths in plot_data.values():
            sample_lengths = set(lengths.keys())
            if len(sample_lengths) > 1:
                multiple_lengths = True
            unique_lengths.update(sample_lengths)
        
        if not multiple_lengths:
            if len(unique_lengths) == 1:
                msg = self.all_same_html
            else:
                msg = self.all_same_within_samples_html
            html = msg.format(",".join(unique_lengths))
        else:
            plot_config = {
                'id': 'fastqc_sequence_length_distribution_plot',
                'title': 'Sequence Length Distribution',
                'ylab': 'Read Count',
                'xlab': 'Sequence Length (bp)',
                'ymin': 0,
                'yMinTickInterval': 0.1,
                'xDecimals': False,
                'colors': self.get_status_cols('sequence_length_distribution'),
                'tt_label': '<b>{point.x} bp</b>: {point.y}',
            }
            
            html = self.plot_html.format(
                linegraph.plot(plot_data, plot_config))
        
        return html



        # Prep the data
        data = dict()
        for s_name in self.fastqc_data:
            bs = self.fastqc_data[s_name]['basic_statistics']
            data[s_name] = {
                'percent_gc': bs['%GC'],
                'avg_sequence_length': bs['avg_sequence_length'],
                'total_sequences': bs['Total Sequences'],
            }
            try:
                data[s_name]['percent_duplicates'] = 100 - bs['total_deduplicated_percentage']
            except KeyError:
                pass # Older versions of FastQC don't have this
            # Add count of fail statuses
            num_statuses = 0
            num_fails = 0
            for s in self.fastqc_data[s_name]['statuses'].values():
                num_statuses += 1
                if s == 'fail':
                    num_fails += 1
            data[s_name]['percent_fails'] = (float(num_fails)/float(num_statuses))*100.0

        # Are sequence lengths interesting?
        seq_lengths = [x['avg_sequence_length'] for x in data.values()]
        hide_seq_length = False if max(seq_lengths) - min(seq_lengths) > 10 else True

        
        headers['percent_duplicates'] = {
            'title': '% Dups',
            'description': '% Duplicate Reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlGn-rev',
            'format': '{:.1f}%'
        }
        headers['percent_gc'] = {
            'title': '% GC',
            'description': 'Average % GC Content',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'Set1',
            'format': '{:.0f}%'
        }
        headers['avg_sequence_length'] = {
            'title': 'Length',
            'description': 'Average Sequence Length (bp)',
            'min': 0,
            'suffix': 'bp',
            'scale': 'RdYlGn',
            'format': '{:.0f}',
            'hidden': hide_seq_length
        }
        headers['percent_fails'] = {
            'title': '% Failed',
            'description': 'Percentage of modules failed in FastQC report (includes those not plotted here)',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'Reds',
            'format': '{:.0f}%',
            'hidden': True
        }
        
        self.general_stats_addcols(data, headers)
    
    
    
    
    def parse_cutadapt_logs(self, f):
        """ Go through log file looking for cutadapt output """
        fh = f['f']
        regexes = {
            '1.7': {
                'bp_processed': "Total basepairs processed:\s*([\d,]+) bp",
                'bp_written': "Total written \(filtered\):\s*([\d,]+) bp",
                'quality_trimmed': "Quality-trimmed:\s*([\d,]+) bp",
                'r_processed': "Total reads processed:\s*([\d,]+)",
                'r_with_adapters': "Reads with adapters:\s*([\d,]+)"
            },
            '1.6': {
                'r_processed': "Processed reads:\s*([\d,]+)",
                'bp_processed': "Processed bases:\s*([\d,]+) bp",
                'r_trimmed': "Trimmed reads:\s*([\d,]+)",
                'quality_trimmed': "Quality-trimmed:\s*([\d,]+) bp",
                'bp_trimmed': "Trimmed bases:\s*([\d,]+) bp",
                'too_short': "Too short reads:\s*([\d,]+)",
                'too_long': "Too long reads:\s*([\d,]+)",
            }
        }
        s_name = None
        cutadapt_version = '1.7'
        for l in fh:
            # New log starting
            if 'cutadapt' in l:
                s_name = None
                c_version = re.match(r'This is cutadapt ([\d\.]+)', l)
                if c_version:
                    try:
                        assert(StrictVersion(c_version.group(1)) <= StrictVersion('1.6'))
                        cutadapt_version = '1.6'
                    except:
                        cutadapt_version = '1.7'
                c_version_old = re.match(r'cutadapt version ([\d\.]+)', l)
                if c_version_old:
                    try:
                        assert(StrictVersion(c_version.group(1)) <= StrictVersion('1.6'))
                        cutadapt_version = '1.6'
                    except:
                        # I think the pattern "cutadapt version XX" is only pre-1.6?
                        cutadapt_version = '1.6'
            # Get sample name from end of command line params
            if l.startswith('Command line parameters'):
                s_name = l.split()[-1]
                # Manage case where sample name is '-' (reading from stdin)
                if s_name == '-':
                    s_name = f['s_name']
                else:
                    s_name = self.clean_s_name(s_name, f['root'])
                if s_name in self.cutadapt_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                self.cutadapt_data[s_name] = dict()
                self.cutadapt_length_counts[s_name] = dict()
                self.cutadapt_length_exp[s_name] = dict()
                self.cutadapt_length_obsexp[s_name] = dict()

            if s_name is not None:
                self.add_data_source(f, s_name)

                # Search regexes for overview stats
                for k, r in regexes[cutadapt_version].items():
                    match = re.search(r, l)
                    if match:
                        self.cutadapt_data[s_name][k] = int(match.group(1).replace(',', ''))

                # Histogram showing lengths trimmed
                if 'length' in l and 'count' in l and 'expect' in l:
                    # Nested loop to read this section while the regex matches
                    for l in fh:
                        r_seqs = re.search("^(\d+)\s+(\d+)\s+([\d\.]+)", l)
                        if r_seqs:
                            a_len = int(r_seqs.group(1))
                            self.cutadapt_length_counts[s_name][a_len] = int(r_seqs.group(2))
                            self.cutadapt_length_exp[s_name][a_len] = float(r_seqs.group(3))
                            if float(r_seqs.group(3)) > 0:
                                self.cutadapt_length_obsexp[s_name][a_len] = float(r_seqs.group(2)) / float(r_seqs.group(3))
                            else:
                                # Cheating, I know. Infinity is difficult to plot.
                                self.cutadapt_length_obsexp[s_name][a_len] = float(r_seqs.group(2))
                        else:
                            break

        # Calculate a few extra numbers of our own
        for s_name, d in self.cutadapt_data.items():
            if 'bp_processed' in d and 'bp_written' in d:
                self.cutadapt_data[s_name]['percent_trimmed'] = (float(d['bp_processed'] - d['bp_written']) / d['bp_processed']) * 100
            elif 'bp_processed' in d and 'bp_trimmed' in d:
                self.cutadapt_data[s_name]['percent_trimmed'] = ((float(d.get('bp_trimmed', 0)) + float(d.get('quality_trimmed', 0))) / d['bp_processed']) * 100



    









    def parse_fastqc_report(self, file_contents, s_name=None, f=None):
        """ Takes contents from a fastq_data.txt file and parses out required
        statistics and data. Returns a dict with keys 'stats' and 'data'.
        Data is for plotting graphs, stats are for top table. """

        # Make the sample name from the input filename if we find it
        fn_search = re.search(r"Filename\s+(.+)", file_contents)
        if fn_search:
            s_name = self.clean_s_name(fn_search.group(1) , f['root'])

        if s_name in self.fastqc_data.keys():
            log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
        self.add_data_source(f, s_name)
        self.fastqc_data[s_name] = { 'statuses': dict() }

        # Parse the report
        section = None
        s_headers = None
        self.dup_keys = []
        for l in file_contents.splitlines():
            if l == '>>END_MODULE':
                section = None
                s_headers = None
            elif l.startswith('>>'):
                (section, status) = l[2:].split("\t", 1)
                section = section.lower().replace(' ', '_')
                self.fastqc_data[s_name]['statuses'][section] = status
            elif section is not None:
                if l.startswith('#'):
                    s_headers = l[1:].split("\t")
                    # Special case: Total Deduplicated Percentage header line
                    if s_headers[0] == 'Total Deduplicated Percentage':
                        self.fastqc_data[s_name]['basic_statistics'].append({
                            'measure': 'total_deduplicated_percentage',
                            'value': float(s_headers[1])
                        })
                    else:
                        # Special case: Rename dedup header in old versions of FastQC (v10)
                        if s_headers[1] == 'Relative count':
                            s_headers[1] = 'Percentage of total'
                        s_headers = [s.lower().replace(' ', '_') for s in s_headers]
                        self.fastqc_data[s_name][section] = list()

                elif s_headers is not None:
                    s = l.split("\t")
                    row = dict()
                    for (i, v) in enumerate(s):
                        v.replace('NaN','0')
                        try:
                            v = float(v)
                        except ValueError:
                            pass
                        row[s_headers[i]] = v
                    self.fastqc_data[s_name][section].append(row)
                    # Special case - need to remember order of duplication keys
                    if section == 'sequence_duplication_levels':
                        try:
                            self.dup_keys.append(float(s[0]))
                        except ValueError:
                            self.dup_keys.append(s[0])

        # Tidy up the Basic Stats
        self.fastqc_data[s_name]['basic_statistics'] = {d['measure']: d['value'] for d in self.fastqc_data[s_name]['basic_statistics']}

        # Calculate the average sequence length (Basic Statistics gives a range)
        length_bp = 0
        total_count = 0
        for d in self.fastqc_data[s_name].get('sequence_length_distribution', {}):
            length_bp += d['count'] * self.avg_bp_from_range(d['length'])
            total_count += d['count']
        if total_count > 0:
            self.fastqc_data[s_name]['basic_statistics']['avg_sequence_length'] = length_bp / total_count

    # These stats are not yet implemented in Atropos.
    #
    # def seq_dup_levels_plot (self):
    #     """ Create the HTML for the Sequence Duplication Levels plot """
    #
    #     data = dict()
    #     for s_name in self.fastqc_data:
    #         try:
    #             d = {d['duplication_level']: d['percentage_of_total'] for d in self.fastqc_data[s_name]['sequence_duplication_levels']}
    #             data[s_name] = OrderedDict()
    #             for k in self.dup_keys:
    #                 try:
    #                     data[s_name][k] = d[k]
    #                 except KeyError:
    #                     pass
    #         except KeyError:
    #             pass
    #     if len(data) == 0:
    #         log.debug('sequence_length_distribution not found in FastQC reports')
    #         return None
    #
    #     pconfig = {
    #         'id': 'fastqc_sequence_duplication_levels_plot',
    #         'title': 'Sequence Duplication Levels',
    #         'categories': True,
    #         'ylab': '% of Library',
    #         'xlab': 'Sequence Duplication Level',
    #         'ymax': 100,
    #         'ymin': 0,
    #         'yMinTickInterval': 0.1,
    #         'colors': self.get_status_cols('sequence_duplication_levels'),
    #         'tt_label': '<b>{point.x}</b>: {point.y:.1f}%',
    #     }
    #
    #     self.sections.append({
    #         'name': 'Sequence Duplication Levels',
    #         'anchor': 'fastqc_sequence_duplication_levels',
    #         'content': '<p>The relative level of duplication found for every sequence. ' +
    #                     'See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html" target="_bkank">FastQC help</a>.</p>' +
    #                     linegraph.plot(data, pconfig)
    #     })
    #
    # def overrepresented_sequences (self):
    #     """Sum the percentages of overrepresented sequences and display them in a bar plot"""
    #
    #     data = dict()
    #     for s_name in self.fastqc_data:
    #         data[s_name] = dict()
    #         try:
    #             max_pcnt   = max( [ float(d['percentage']) for d in self.fastqc_data[s_name]['overrepresented_sequences']] )
    #             total_pcnt = sum( [ float(d['percentage']) for d in self.fastqc_data[s_name]['overrepresented_sequences']] )
    #             data[s_name]['total_overrepresented'] = total_pcnt
    #             data[s_name]['top_overrepresented'] = max_pcnt
    #             data[s_name]['remaining_overrepresented'] = total_pcnt - max_pcnt
    #         except KeyError:
    #             if self.fastqc_data[s_name]['statuses']['overrepresented_sequences'] == 'pass':
    #                 data[s_name]['total_overrepresented'] = 0
    #                 data[s_name]['top_overrepresented'] = 0
    #                 data[s_name]['remaining_overrepresented'] = 0
    #             else:
    #                 log.debug("Couldn't find data for {}, invalid Key".format(s_name))
    #
    #     cats = OrderedDict()
    #     cats['top_overrepresented'] = { 'name': 'Top over-represented sequence' }
    #     cats['remaining_overrepresented'] = { 'name': 'Sum of remaining over-represented sequences' }
    #
    #     # Config for the plot
    #     pconfig = {
    #         'id': 'fastqc_overrepresented_sequencesi_plot',
    #         'title': 'Overrepresented sequences',
    #         'ymin': 0,
    #         'yCeiling': 100,
    #         'yMinRange': 20,
    #         'tt_decimals': 2,
    #         'tt_suffix': '%',
    #         'tt_percentages': False,
    #         'ylab_format': '{value}%',
    #         'cpswitch': False,
    #         'ylab': 'Percentage of Total Sequences'
    #     }
    #
    #     # Check if any samples have more than 1% overrepresented sequences, else don't make plot.
    #     if max([ x['total_overrepresented'] for x in data.values()]) < 1:
    #         plot_html = '<div class="alert alert-info">{} samples had less than 1% of reads made up of overrepresented sequences</div>'.format(len(data))
    #     else:
    #         plot_html = bargraph.plot(data, cats, pconfig)
    #
    #     self.sections.append({
    #         'name': 'Overrepresented sequences',
    #         'anchor': 'fastqc_overrepresented_sequences',
    #         'content': '<p> The total amount of overrepresented sequences found in each library. ' +
    #                 'See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/9%20Overrepresented%20Sequences.html" target="_bkank">FastQC help for further information</a>.</p>'
    #                 + plot_html
    #         })
    #
    # def adapter_content_plot (self):
    #     """ Create the HTML for the FastQC adapter plot """
    #
    #     data = dict()
    #     for s_name in self.fastqc_data:
    #         try:
    #             for d in self.fastqc_data[s_name]['adapter_content']:
    #                 pos = self.avg_bp_from_range(d['position'])
    #                 for r in self.fastqc_data[s_name]['adapter_content']:
    #                     pos = self.avg_bp_from_range(r['position'])
    #                     for a in r.keys():
    #                         k = "{} - {}".format(s_name, a)
    #                         if a != 'position':
    #                             try:
    #                                 data[k][pos] = r[a]
    #                             except KeyError:
    #                                 data[k] = {pos: r[a]}
    #         except KeyError:
    #             pass
    #     if len(data) == 0:
    #         log.debug('adapter_content not found in FastQC reports')
    #         return None
    #
    #     # Lots of these datasets will be all zeros.
    #     # Only take datasets with > 0.1% adapter contamination
    #     data = {k:d for k,d in data.items() if max(data[k].values()) >= 0.1 }
    #
    #     pconfig = {
    #         'id': 'fastqc_adapter_content_plot',
    #         'title': 'Adapter Content',
    #         'ylab': '% of Sequences',
    #         'xlab': 'Position',
    #         'yCeiling': 100,
    #         'yMinRange': 5,
    #         'ymin': 0,
    #         'xDecimals': False,
    #         'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}%',
    #         'hide_empty': True,
    #         'yPlotBands': [
    #             {'from': 20, 'to': 100, 'color': '#e6c3c3'},
    #             {'from': 5, 'to': 20, 'color': '#e6dcc3'},
    #             {'from': 0, 'to': 5, 'color': '#c3e6c3'},
    #         ],
    #     }
    #
    #     if len(data) > 0:
    #         plot_html = linegraph.plot(data, pconfig)
    #     else:
    #         plot_html = '<div class="alert alert-warning">No samples found with any adapter contamination > 0.1%</div>'
    #
    #     # Note - colours are messy as we've added adapter names here. Not
    #     # possible to break down pass / warn / fail for each adapter, which
    #     # could lead to misleading labelling (fails on adapter types with
    #     # little or no contamination)
    #
    #     self.sections.append({
    #         'name': 'Adapter Content',
    #         'anchor': 'fastqc_adapter_content',
    #         'content': '<p>The cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position. ' +
    #                     'See the <a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html" target="_bkank">FastQC help</a>. ' +
    #                     'Only samples with &ge; 0.1% adapter contamination are shown.</p>' + plot_html
    #     })
