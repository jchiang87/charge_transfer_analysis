import matplotlib
import lsst.charge_transfer_analysis as cta

matplotlib.pyplot.ion()

datapath = '/nfs/farm/g/lsst/u1/jobHarness/jh_archive/ITL-CCD/ITL-3800C-065/cte_offline/v0/4557'
sensor_id = 'ITL-3800C-065'

results = cta.fit_superflats(datapath, sensor_id)
results.plots.savefig('%s_overscan_fits.png' % sensor_id)
results.fit_results.to_pickle('%s_fit_results.pkl' % sensor_id)
