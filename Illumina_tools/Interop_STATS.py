# Run in 
# conda activate Illumina_interop

## Import runfolder
import argparse

parser = argparse.ArgumentParser(description='Process Illumina raw data.')
parser.add_argument('runfolder', action='store', help='full path to run directory')
args = parser.parse_args()
runfolder = args.runfolder
print(">> Analyzing Runfolder:")
print(runfolder)

# Import libs
from interop import py_interop_run_metrics, py_interop_run, py_interop_summary

## Run Info:
run_info = py_interop_run.info()
run_info.read(runfolder)
print(">> Run Info: ")
print("- Number of cycles %d" % run_info.total_cycles())
print("- Is paired end: %s" % run_info.is_paired_end())

            
run_metrics = py_interop_run_metrics.run_metrics()
valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
py_interop_run_metrics.list_summary_metrics_to_load(valid_to_load)

run_folder = run_metrics.read(runfolder, valid_to_load)
summary = py_interop_summary.run_summary()
py_interop_summary.summarize_run_metrics(run_metrics, summary)


## Yield
import pandas as pd

columns = ( ('Yield Total (G bases)', 'yield_g'), ('% Aligned PhiX', 'percent_aligned'), ('Percent Q30', 'percent_gt_q30') )
rows = [('Non-Indexed Total', summary.nonindex_summary()), ('Total', summary.total_summary())]
d = []

for label, func in columns:
    d.append( (label, pd.Series([getattr(r[1], func)() for r in rows], index=[r[0] for r in rows])))

df = pd.DataFrame.from_dict(dict(d))
print("\n>> Yield:")
print(df)
df.round(2).to_csv(str(runfolder) + '/Interop_Stats_Yield.csv')


## Yield per read
columns = ( ('Yield Total (G bases)', 'yield_g'),  ('% Aligned PhiX', 'percent_aligned'))
rows = [("Read %s%d"%("(I)" if summary.at(i).read().is_index()  else " ", summary.at(i).read().number()), summary.at(i).summary()) for i in range(summary.size())]
d = []
for label, func in columns:
    d.append( (label, pd.Series([getattr(r[1], func)() for r in rows], index=[r[0] for r in rows])))

df2 = pd.DataFrame.from_dict(dict(d))
print("\n>> Yield pr read:")
print(df2)
df2.round(2).to_csv(str(runfolder) + '/Interop_Stats_Yield.per.read.csv')

## Density
def format_value(val):
    if hasattr(val, 'mean'):
        return val.mean()
    else:
        return val

columns = ( ('Lane', 'lane'), ('Tiles', 'tile_count'), ('Density (K/mm2)', 'density'), ('Total Clusters', 'reads'), ('Clusters PF', 'reads_pf'))

read = 0

rows = [summary.at(read).at(lane) for lane in range(summary.lane_count())]
d = []

for label, func in columns:
    d.append( (label, pd.Series([format_value(getattr(r, func)()) for r in rows])))


df = pd.DataFrame.from_dict(dict(d))
print("\n>> Lane Stats:")
print(df)
df.round(2).to_csv(str(runfolder) + '/Interop_Stats_Lanes.csv')
