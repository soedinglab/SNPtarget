from setuptools import setup, find_packages

__name__ = "enh_scripts"
__version__ = "0.2.0"
__description__ = ("helper scripts for the prediction of enhancer-gene "
                   "interactions")
__author__ = "Christian Roth"
__author_email__ = "christian.roth@mpibpc.mpg.de"

install_deps = [
    "pandas",
    "ngsbiotools",
    "pyivtree",
]
# we don't want to risk updating certain scientific packages
try:
    import numpy as np
    np.__version__
except:
    install_deps.append("numpy")

try:
    import scipy as sp
    sp.__version__
except:
    install_deps.append("scipy")

setup(
    name=__name__,
    version=__version__,
    description=__description__,
    author=__author__,
    author_email=__author_email__,
    packages=find_packages(),
    install_requires=install_deps,
    test_suite='nose.collector',
    entry_points={
        'console_scripts': [
            'gc2trscr = enh_scripts.gc2trscr:main',
            'trscr2reg = enh_scripts.trscr2reg:main',
            'sheff_bgsampling = enh_scripts.sheff_bgsampling:main',
            'calc_dist = enh_scripts.calc_dist:main',
            'calc_posterior = enh_scripts.calc_posterior:main',
            'posterior_pval = enh_scripts.posterior_pval:main',
            'noise_control = enh_scripts.noise_control:main',
            'split_dhs = enh_scripts.split_dhs:main',
            'split_expr = enh_scripts.split_expr:main',
            'rao_hic_to_matrix = enh_scripts.rao_hic_to_matrix:main',
            'extract_access = enh_scripts.extract_access:main',
            'extract_cuts_general = enh_scripts.extract_cuts_general:main',
            'extract_cuts_intervals = enh_scripts.extract_cuts_intervals:main',
            'extract_hic = enh_scripts.extract_hic:main',
            'hic_pval = enh_scripts.hic_pval:main',
            'dbsnp2loc = enh_scripts.dbsnp2loc:main',
            'extract_cuts = enh_scripts.extract_cuts:main',
            'extract_tss_cuts = enh_scripts.extract_tss_cuts:main',
            'normalize_cuts = enh_scripts.normalize_cuts:main',
            'tgvcf2snpdb = enh_scripts.tgvcf2snpdb:main',
            'snp_access_mask = enh_scripts.snp_access_mask:main',
            'corr_pval = enh_scripts.corr_pval:main',
            'corr_pval_test = enh_scripts.corr_pval_test:main',
            'expr_error_estim = enh_scripts.expr_error_estim:main',
            'plot_pval_score = enh_scripts.plot_pval_score:main',
            'fit_pval_bf = enh_scripts.fit_pval_bf:main',
            'fit_dist_bf = enh_scripts.fit_dist_bf:main',
            'fit_hic_bf = enh_scripts.fit_hic_bf:main',
            'estimate_eta0 = enh_scripts.estimate_eta0:main',
            'skip_counter = enh_scripts.skip_counter:main',
            'list_eg_pairs = enh_scripts.list_enh_gene_pairs:main',
            'subsample_dist_bf = enh_scripts.subsample_dist_bf:main',
            'filter_promoters = enh_scripts.filter_promoters:main',
            'select_tss = enh_scripts.select_tss:main',
            'extract_snpids = enh_scripts.extract_snpid_by_cellline:main',
            'annotate_with_ia = enh_scripts.annotate_with_ia:main',
            'filter_by_snpids = enh_scripts.filter_by_snpids:main',
            'ld_snp_finder = enh_scripts.saikat_snp_finder:main',
        ]
    },
)
