"""MultiQC module that reads JSON output from Atropos.

Note: this code is mostly borrowed from the Cutadapt [1] and FastQC [2] MultiQC
modules.
1. https://github.com/ewels/MultiQC/blob/master/multiqc/modules/fastqc/fastqc.py
2. https://github.com/ewels/MultiQC/blob/master/multiqc/modules/cutadapt/cutadapt.py
"""
from __future__ import print_function, division, absolute_import
from collections import OrderedDict, defaultdict
import logging
import io
import math
from numpy import mean, multiply, divide, cumsum
import operator
import os
import json
import re

from multiqc import config
from multiqc.plots import linegraph
from multiqc.modules.base_module import BaseMultiqcModule

## Hard-coded URLs
# TODO: eventually these need pointers to specific sections of the docs

ATROPOS_GITHUB_URL = "https://github.com/jdidion/atropos"
ATROPOS_DOC_URL = "http://atropos.readthedocs.org/en/latest/guide.html"

## Assets

ATROPOS_CSS = {
    'assets/css/multiqc_atropos.css' :
    os.path.join(os.path.dirname(__file__), 'assets', 'css', 'multiqc_atropos.css')
}
ATROPOS_JS = {
    'assets/js/multiqc_atropos.js' :
    os.path.join(os.path.dirname(__file__), 'assets', 'js', 'multiqc_atropos.js')
}

class Status(object):
    def __init__(self, name, val, color):
        self.name = name
        self.val = val
        self.color = color
    
    def __lt__(self, other):
        return self.val < other.val
    
    def __eq__(self, other):
        return self.val == other.val
    
    def __repr__(self):
        return self.name

FAIL = Status('fail', 0, '#5cb85c')
WARN = Status('warn', 1, '#f0ad4e')
PASS = Status('pass', 2, '#d9534f')

DEFAULT_COLOR = '#999'

ATROPOS_TRIMMED_LENGTH = """
<p>This plot shows the number of reads with certain lengths of adapter trimmed.
Obs/Exp shows the raw counts divided by the number expected due to sequencing
errors. A defined peak may be related to adapter length. See the
<a href="{}" target="_blank">Atropos documentation</a> for more information on
how these numbers are generated.</p>""".format(ATROPOS_DOC_URL)

ATROPOS_PASSFAILS = '<script type="text/javascript">atropos_passfails = {};</script>'

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
        self.atropos_general_data = OrderedDict()
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
        # Atropos summary files are JSON files (with .json extension) and
        # always have '"program": "Atropos"' as the first key-value pair.
        if 'atropos' in config.sp:
            patterns = config.sp['atropos']
        else:
            patterns = dict(fn='*.json', contents='"program": "Atropos"')
        for file_dict in self.find_log_files(patterns, filehandles=True):
            fileobj = file_dict['f']
            data = json.load(fileobj)
            self.add_atropos_data(data, fileobj, file_dict)
        
        num_samples = len(self.atropos_sample_ids)
        if num_samples == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning
        else:
            log.info("Found {} reports".format(num_samples))
    
    def add_atropos_data(self, data=None, data_source=None, file_dict=None):
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
        if file_dict:
            self.add_data_source(file_dict, sample_id)
        
        self.atropos_general_data[sample_id] = data['derived'].copy()
        self.atropos_general_data[sample_id].update(dict(
            (key, data[key])
            for key in (
                'total_record_count',
                'total_bp_counts')))
        
        if 'trim' in data:
            self.atropos_trim_data[sample_id] = data['trim']
            self.atropos_general_data[sample_id]['fraction_bp_trimmed'] = \
                1.0 - data['trim']['formatters']['fraction_total_bp_written']
            for name, mod_dict in data['trim']['modifiers'].items():
                if name in ('AdapterCutter', 'InsertAdapterCutter'):
                    self.atropos_general_data[sample_id]['fraction_records_with_adapters'] = \
                        mod_dict['fraction_records_with_adapters']
                    break
        
        pairing = data['input']['input_read']
        
        def add_stats_data(phase, source, data):
            self.atropos_qc_data[phase][sample_id] = OrderedDict()
            source_files = source.split(',')
            if pairing == 3:
                self.atropos_qc_data[phase][sample_id] = [None, None]
                for idx in range(2):
                    read = 'read{}'.format(idx+1)
                    if read in data:
                        self.atropos_qc_data[phase][sample_id][idx] = data[read]
                        self.atropos_qc_data[phase][sample_id][idx]['source'] = source_files[idx]
            else:
                read = 'read{}'.format(pairing)
                data[read]['source'] = source_files[pairing-1]
                self.atropos_qc_data[phase][sample_id] = (data[read], None)
        
        if 'pre' in data:
            for source, pre_data in data['pre'].items():
                add_stats_data('pre', source, pre_data)
        
        # TODO: handle multiplexed output
        if 'post' in data and 'NoFilter' in data['post']:
            for source, post_data in data['post']['NoFilter'].items():
                add_stats_data('post', source, post_data)
    
    def atropos_report(self):
        if not self.atropos_sample_ids:
            log.debug("No reports to process; atropos_report raising UserWarning")
            raise UserWarning
        
        # Add to general stats
        if self.atropos_general_data:
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
        headers['input_format'] = {
            'title': 'Format',
            'description': 'Input File Format'
        }
        headers['total_record_count'] = {
            'title': 'M Seqs',
            'description': 'Total Sequences (millions)',
            'min': 0,
            'scale': 'Blues',
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }
        headers['total_bp_counts'] = {
            'title': 'M bp',
            'description': 'Total Base Pairs (millions)',
            'min': 0,
            'scale': 'Blues',
            'modify': lambda x: x / 1000000,
            'shared_key': 'base_count'
        }
        headers['avg_sequence_length'] = {
            'title': 'Length',
            'description': 'Average Sequence Length (bp)',
            'min': 0,
            'suffix': 'bp',
            'scale': 'RdYlGn',
            'format': '{:.0f}'
        }
        headers['percent_fails_pre'] = {
            'title': '% Failed (pre-trimming)',
            'description': 'Percentage of Atropos QC modules failed (pre-trimming)',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'Reds',
            'format': '{:.0f}%',
            'hidden': True
        }
        headers['percent_fails_post'] = {
            'title': '% Failed (post-trimming)',
            'description': 'Percentage of Atropos QC modules failed (post-trimming)',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'Reds',
            'format': '{:.0f}%',
            'hidden': True
        }
        self.general_stats_addcols(self.atropos_general_data, headers)
    
    def atropos_qc(self, phase):
        # Add statuses to intro. Note that this is slightly different than
        # FastQC: Atropos reports the relevant statistic, and the threshold
        # for pass/warn/fail is configured in MutliQC (defaulting to
        # the thresholds defined in FastQC).
        statuses = {}
        fails = defaultdict(int)
        phase_data = self.atropos_qc_data[phase]
        def get_section_data(section_name):
            section_data = {}
            for sample_id, (read1_data, read2_data) in phase_data.items():
                sample_data = [None, None]
                if read1_data:
                    sample_data[0] = read1_data[section_name]
                if read2_data:
                    sample_data[1] = read2_data[section_name]
                section_data[sample_id] = sample_data
            return section_data
        for section in self.sections:
            section_name = "{}_{}".format(phase, section.name)
            section_data = get_section_data(section.name)
            status, html = section(phase, section_data)
            statuses[section_name] = status
            for sample_id, status in statuses[section_name].items():
                if status == FAIL:
                    fails[sample_id] += 1
            self.section_html.append(html)
        
        for sample_id, num_fails in fails.items():
            self.atropos_general_data[sample_id]['percent_fails_{}'.format(phase)] = \
                num_fails * 100 / len(fails)
        
        self.intro += ATROPOS_PASSFAILS.format(json.dumps(statuses))
    
    def atropos_trim(self):
        """ Take the parsed stats from the Atropos report and add it to the
        basic stats table at the top of the report.
        """
        headers = {}
        headers['fraction_bp_trimmed'] = {
            'title': '% Trimmed',
            'description': '% Total Base Pairs trimmed',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlBu-rev',
            'format': '{:.1f}%'
        }
        header['fraction_records_with_adapters'] = {
            'title': '% Reads w/ Adapters',
            'description': '% Total Reads with Adapters',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlBu-rev',
            'format': '{:.1f}%'
        }
        self.general_stats_addcols(self.atropos_general_data, headers)

        """ Generate the trimming length plot.
        """
        
        html = ATROPOS_TRIMMED_LENGTH
        
        plot_config = {
            'id': 'atropos_plot',
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
            plot_config)

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
    return d

def weighted_lower_quantile(vals, freqs, quantile):
    for i, frac in enumerate(freqs):
        if frac > quantile:
            return vals[i-1] if i > 0 else vals[0]
    return values[-1]

def weighted_median(vals, counts_cumsum):
    total = counts_cumsum[-1]
    if total == 0:
        return None
    mid1 = mid2 = (total // 2) + 1
    if total % 2 == 0:
        mid1 -= 1
    val1 = val2 = None
    for i, val in enumerate(counts_cumsum):
        if val1 is None and mid1 <= val:
            val1 = vals[i]
        if mid2 <= val:
            val2 = vals[i]
            break
    return float(val1 + val2) / 2

def weighted_mean(vals, counts):
    return sum(multiply(vals, counts)) / sum(counts)

def hist_to_means(values, base_counts):
    return dict(
        (pos, weighted_mean(values, counts))
        for pos, counts in base_counts)

def get_status_cols(statuses):
    """Helper function - returns a list of colours according to the
    status of this module for each sample.
    """
    return ordered_dict(
        (sample_id, status.color)
        for sample_id, status in statuses.items())

## QC sections

class Section(object):
    #headers = {}
    
    def __call__(self, phase, data):
       statuses = self.get_statuses(phase, data)
       #if section.headers:
       #    self.general_stats_addcols(section_data, section.headers)
       html = self.plot(statuses, data)
       return (statuses, html)
    
    def get_statuses(self, phase, data):
        return ordered_dict(
            (sample_id, self.get_status(phase, *sample_data))
            for sample_id, sample_data in data.items())
    
    def get_status(self, phase, sample_data1, sample_data2=None):
        def get_status_stat(data):
            if data is None:
                return None
            try:
                return data[self.threshold_statistic]
            except:
                return self.compute_statistic(data)
        stat1, stat2 = [
            get_status_stat(data)
            for data in (sample_data1, sample_data2)]
        return self.get_status_for_pair(phase, stat1, stat2)
    
    def compute_statistic(self, data):
        raise NotImplementedError()
    
    def get_status_for_pair(self, phase, stat1, stat2=None):
        status1 = self.get_status_for(phase, stat1)
        if stat2:
            status2 = self.get_status_for(phase, stat2)
            return min(status1, status2)
        else:
            return status1
    
    def get_status_for(self, phase, stat):
        # TODO: Add configuration for thresholds.
        for status, threshold in zip((FAIL, WARN), self.default_threshold):
            if self.compare(stat, threshold):
                return status
        else:
            return PASS
    
    def plot(self, statuses, data):
        return {
            'name': self.display,
            'anchor': self.anchor,
            'content': self.get_html(statuses, data)
        }
    
    def get_html(self, statuses, data):
        plot_data1, plot_data2, extra = self.get_plot_data(data)
        return self.html.format(*self.get_html_variables(
            statuses, plot_data1, plot_data2, **extra))
    
    def get_plot_data(self, data):
        plot_data1 = OrderedDict()
        plot_data2 = OrderedDict()
        for sample_id, (data1, data2) in data.items():
            plot_data1[sample_id] = self.get_sample_plot_data(data1)
            if data2:
                plot_data2[sample_id] = self.get_sample_plot_data(data2)
        return (plot_data1, plot_data2, {})
    
    def get_sample_plot_data(self, data):
        raise NotImplementedError()
    
    def get_html_variables(self, statuses, plot_data1, plot_data2, **extra):
        content1 = self.get_plot(statuses, plot_data1, **extra)
        if plot_data2:
            content2 = self.get_plot(statuses, plot_data2, **extra)
            return (self.wrap_plot_pair(content1, content2),)
        else:
            return (content1,)
    
    def get_plot(self, statuses, plot_data):
        raise self.plot_type.plot(
            plot_data, self.get_plot_config(statuses))
    
    def wrap_plot_pair(self, content1, content2):
        return '<div class="pair-plot-wrapper">' + content1 + content2 + '<\div>'

## Static HTML/templates

class PerBaseQuality(Section):
    name = 'base_qualities'
    display = 'Sequence Quality Histograms'
    anchor = 'atropos_per_base_sequence_quality'
    html = """
<p>The mean quality value across each base position in the read. See the
<a href="{}" target="_bkank">Atropos help</a>.</p>
{{}}""".format(ATROPOS_DOC_URL)
    # TODO: Add boxplots as in FastQC output.
    plot_type = linegraph
    
    def get_plot_config(self, statuses):
        return {
            'id': 'atropos_per_base_sequence_quality_plot',
            'title': 'Mean Quality Scores',
            'ylab': 'Phred Score',
            'xlab': 'Position (bp)',
            'ymin': 0,
            'xDecimals': False,
            'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}',
            'colors': get_status_cols(statuses),
            'yPlotBands': [
                {'from': 28, 'to': 100, 'color': '#c3e6c3'},
                {'from': 20, 'to': 28, 'color': '#e6dcc3'},
                {'from': 0, 'to': 20, 'color': '#e6c3c3'},
            ]
        }
    
    def compute_statistic(self, data):
        # base_qualities is in the form
        # [[qualities], [[base1_counts], [base2_counts], ...]]
        quals = data[0]
        min_lower_quartile = None
        min_median = None
        for pos, base_counts in data[1]:
            counts_cumsum = cumsum(base_counts)
            lower_quartile = weighted_lower_quantile(
                quals,
                (float(c) / counts_cumsum[-1] for c in base_counts),
                0.25)
            if min_lower_quartile is None or lower_quartile < min_lower_quartile:
                min_lower_quartile = lower_quartile
            median = weighted_median(quals, counts_cumsum)
            if min_median is None or median < min_median:
                min_median = median
        return (min_lower_quartile, min_median)
    
    def get_status_for(self, phase, stat):
        min_lower_quartile, min_median = stat
        if min_lower_quartile < 5 or min_median < 20:
            return FAIL
        if min_lower_quartile < 10 or min_median < 25:
            return WARN
        return PASS
    
    def get_sample_plot_data(self, data):
        return hist_to_means(*data)

class PerTileQuality(Section):
    name = 'tile_sequence_quality'
    compare = operator.lt
    threshold_statistic = 'min_tile_versus_mean'
    default_thresholds = (-5, -2)
    # TODO

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
    
    def get_plot_config(self, statuses):
        return {
            'id': 'atropos_per_sequence_quality_scores_plot',
            'title': 'Per Sequence Quality Scores',
            'ylab': 'Count',
            'xlab': 'Mean Sequence Quality (Phred Score)',
            'ymin': 0,
            'xmin': 0,
            'xDecimals': False,
            'colors': get_status_cols(statuses),
            'tt_label': '<b>Phred {point.x}</b>: {point.y} reads',
            'xPlotBands': [
                {'from': 28, 'to': 100, 'color': '#c3e6c3'},
                {'from': 20, 'to': 28, 'color': '#e6dcc3'},
                {'from': 0, 'to': 20, 'color': '#e6c3c3'},
            ]
        }
    
    def get_plot_data(self, data):
        return hist_to_means(data)

class PerBaseContent(Section):
    name = 'bases'
    display = 'Per Base Sequence Content'
    anchor = 'atropos_per_base_sequence_content'
    header_html = """
<p>The proportion of each base position for which each of the four normal DNA
bases has been called. See the <a href="{}" target="_bkank">Atropos help</a>.</p>
<p class="text-primary"><span class="glyphicon glyphicon-info-sign"></span>
Click a heatmap row to see a line plot for that dataset.</p>""".format(ATROPOS_DOC_URL)
    plot_html = """
<div id="atropos_per_base_sequence_content_plot">
    <h5><span class="s_name"><em class="text-muted">rollover for sample name</em></span></h5>
    <div class="atropos_seq_heatmap_key">
        Position: <span id="atropos_seq_heatmap_key_pos">-</span>
        <div><span id="atropos_seq_heatmap_key_t"> %T: <span>-</span></span></div>
        <div><span id="atropos_seq_heatmap_key_c"> %C: <span>-</span></span></div>
        <div><span id="atropos_seq_heatmap_key_a"> %A: <span>-</span></span></div>
        <div><span id="atropos_seq_heatmap_key_g"> %G: <span>-</span></span></div>
    </div>
    <div id="atropos_seq_heatmap_div" class="atropos-overlay-plot">
        <div id="atropos_seq" class="hc-plot">
            <canvas id="atropos_seq_heatmap" height="100%" width="800px" style="width:100%;"></canvas>
        </div>
    </div>
    <div class="clearfix"></div>
</div>
<script type="text/javascript">
    atropos_seq_content_data = {};
    $(function () {{ atropos_seq_content_heatmap(); }});
</script>"""
    
    def get_html(self, statuses, data):
        plot_data1, plot_data2 = self.get_plot_data(data)
        return self.header_html + self.get_html_variables(
            statuses, plot_data1, plot_data2)[0]
    
    def get_sample_plot_data(self, data):
        return dict(data[1])
    
    def get_plot(self, statuses, data):
        """Create the epic HTML for the Atropos sequence content heatmap.
        """
        return self.plot_html.format(json.dumps(data))

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
    
    def get_plot_data(data):
        theoretical_gc = None
        theoretical_gc_name = None
        
        tgc = getattr(config, 'atropos_config', {}).get('atropos_theoretical_gc', None)
        if tgc is not None:
            theoretical_gc_name = os.path.basename(tgc)
            tgc_fn = 'fastqc_theoretical_gc_{}.txt'.format(tgc)
            tgc_path = os.path.join(os.path.dirname(__file__), 'fastqc_theoretical_gc', tgc_fn)
            if not os.path.isfile(tgc_path):
                tgc_path = tgc
            try:
                with io.open (tgc_path, "r", encoding='utf-8') as f:
                    theoretical_gc_raw = f.read()
                    theoretical_gc = []
                    for l in theoretical_gc_raw.splitlines():
                        if '# FastQC theoretical GC content curve:' in l:
                            theoretical_gc_name = l[39:]
                        elif not l.startswith('#'):
                            s = l.split()
                            try:
                                theoretical_gc.append([float(s[0]), float(s[1])])
                            except (TypeError, IndexError):
                                pass
            except IOError:
                log.warn("Couldn't open FastQC Theoretical GC Content file {}".format(tgc_path))
        
        if theoretical_gc is None:
            means = []
            sds = []
            for values in data.values():
                dist = values['sequence_gc']['dist']
                means.append(dist['mean'])
                sds.append(dist['stdev'])
            nd = NormalDistribution(mean(dist['mean']), mean(dist['stdev']))
            theoretical_gc = [nd(i) * total for i in range(101)]
            theoretical_gc_name = "Empirical"
        
        return super().get_html_variables(plot_data1, plot_data2) + (dict(
            theoretical_gc=theoretical_gc,
            theoretical_gc_name=theoretical_gc_name
        ),)
    
    def get_sample_plot_data(self, data):
        def normalize_gc(gc, total):
            return ordered_dict(
                (pct, count * 100 / total)
                for pct, count in gc)
        
        gc = data['sequence_gc']['hist']
        total = sum(gc.values())
        gcnorm = normalize_gc(gc, total)
        return (gc, gcnorm)
    
    def get_html_variables(
            self, statuses, plot_data1, plot_data2, theoretical_gc=None,
            theoretical_gc_name=''):
        return (theoretical_gc_name,) + super().get_html_variables(
            statuses, plot_data1, plot_data2, theoretical_gc=theoretical_gc)
    
    def get_plot(self, statuses, data, theoretical_gc=None):
        plot_data = dict((sample_id, values[0]) for sample_id, values in data.items())
        plot_data_norm = dict((sample_id, values[1]) for sample_id, values in data.items())
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
            'colors': get_status_cols(statuses),
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
        plot_config['extra_series'] = [ [dict(esconfig)], [dict(esconfig)] ]
        plot_config['extra_series'][0][0]['data'] = theoretical_gc
        return linegraph.plot([plot_data_norm, plot_data], plot_config)

class PerBaseN(Section):
    name = 'per_base_N_content'
    display = 'Per Base N Content'
    anchor = 'atropos_per_base_n_content'
    compare = operator.gt
    threshold_statistic = 'frac_N'
    default_thresholds = (0.2, 0.05)
    html = """
<p>The percentage of base calls at each position for which an N was called.
See the <a href="{}" target="_bkank">Atropos help</a>.</p>
{{}}""".format(ATROPOS_DOC_URL)
    plot_type = linegraph
    
    def get_plot_config(self, statuses):
        return {
            'id': 'atropos_per_base_n_content_plot',
            'title': 'Per Base N Content',
            'ylab': 'Percentage N-Count',
            'xlab': 'Position in Read (bp)',
            'yCeiling': 100,
            'yMinRange': 5,
            'ymin': 0,
            'xmin': 0,
            'xDecimals': False,
            'colors': get_status_cols(statuses),
            'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}%',
            'yPlotBands': [
                {'from': 20, 'to': 100, 'color': '#e6c3c3'},
                {'from': 5, 'to': 20, 'color': '#e6dcc3'},
                {'from': 0, 'to': 5, 'color': '#c3e6c3'},
            ]
        }
    
    def get_sample_plot_data(self, data):
        """Create the HTML for the per base N content plot.
        """
        return data['frac_N']

class SequenceLength(Section):
    name = 'sequence_length'
    display = 'Sequence Length Distribution'
    anchor = 'atropos_sequence_length_distribution'
    threshold_statistic = 'length_range'
    html = """
<p>The distribution of fragment sizes (read lengths) found.
See the <a href="{}" target="_bkank">Atropos help</a>.</p>
{{}}""".format(ATROPOS_DOC_URL)
    all_same_html = """
<p>All samples have sequences of a single length ({} bp).</p>"""
    all_same_within_samples_html = """
<p>All samples have sequences of a single length ({} bp).
See the <a href="#general_stats">General Statistics Table</a>.</p>"""
    
    def get_status_for(self, phase, stat):
        """This always returns 'pass' for post-trimming.
        """
        if phase == 'pre':
            min_len, max_len = stat
            if min_len == 0:
                return FAIL
            if min_len != max_len:
                return WARN
        return PASS
    
    def get_html(self, statuses, data):
        def _get_html(plot_data):
            unique_lengths = set()
            multiple_lengths = False
            for hist in data.values():
                if len(hist) > 1:
                    multiple_lengths = True
                unique_lengths.update(h[0] for h in hist)
            if not multiple_lengths:
                if len(unique_lengths) == 1:
                    msg = self.all_same_html
                else:
                    msg = self.all_same_within_samples_html
                return msg.format(",".join(unique_lengths))
            else:
                plot_config = {
                    'id': 'atropos_sequence_length_distribution_plot',
                    'title': 'Sequence Length Distribution',
                    'ylab': 'Read Count',
                    'xlab': 'Sequence Length (bp)',
                    'ymin': 0,
                    'yMinTickInterval': 0.1,
                    'xDecimals': False,
                    'colors': get_status_cols(statuses),
                    'tt_label': '<b>{point.x} bp</b>: {point.y}',
                }
                return self.html.format(linegraph.plot(plot_data, plot_config))
        
        plot_data1, plot_data2, extra = self.get_plot_data(data)
        html1 = _get_html(plot_data1)
        if plot_data2:
            html2 = _get_html(plot_data2)
            return self.wrap_plot_pair(html1, html2)
        else:
            return html1
    
    def get_sample_plot_data(self, data):
        return data['length']

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
    #     plot_config = {
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
    #                     linegraph.plot(data, plot_config)
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
    #     plot_config = {
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
    #         plot_html = bargraph.plot(data, cats, plot_config)
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
    #     plot_config = {
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
    #         plot_html = linegraph.plot(data, plot_config)
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
