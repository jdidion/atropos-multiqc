from pkg_resources import get_distribution
from multiqc.utils import config

__version__ = get_distribution("multiqc_atropos").version
config.multiqc_atropos_version = __version__
