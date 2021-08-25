import pywt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
plt.style.use('seaborn-pastel')
#matplotlib.use('Agg')


import xml.etree.ElementTree as ET
tree = ET.parse('699_PH1_WH4.mzXML')
root = tree.getroot()

import re

basePeakMz1= []
basePeakIntensity1= []
retentionTime1= []
scan1= []

for x in root.findall('.//{http://sashimi.sourceforge.net/schema_revision/mzXML_3.2}scan'):
    x= x.attrib
    scan= int(x['num'])
    scan1.append(scan)
    basePeakMz = x['basePeakMz']
    basePeakMz = re.sub('[^0-9.]+', '', basePeakMz)
    basePeakMz1.append(round(float(basePeakMz), 4))
    basePeakIntensity = x['basePeakIntensity']
    basePeakIntensity = basePeakIntensity.strip()
    basePeakIntensity = re.sub('e0', 'e+0', basePeakIntensity)
    basePeakIntensity1.append(round(float(basePeakIntensity), 4))
    retentionTime = x['retentionTime']
    retentionTime = re.sub('[^0-9.]+', '', retentionTime)
    retentionTime= float(retentionTime)
    retentionTime = round(retentionTime, 4)
    retentionTime1.append(retentionTime)


import pandas as pd
import numpy as np
from scipy.signal import chirp, find_peaks, peak_widths
import scipy.signal

df= pd.DataFrame()
df['scan']= scan1
df['m/z']= basePeakMz1
df['rt']=retentionTime1
df['intensity']= basePeakIntensity1

df.to_csv('raw_peak.csv')


xx= basePeakIntensity1
xx= [abs(x) for x in xx]
x = np.linspace(0, 1, num=2048)
chirp_signal = xx

# def on_plot_hover(event):
#     # Iterating over each data member plotted
#     for curve in plot.get_lines():
#         # Searching which data member corresponds to current mouse position
#         if curve.contains(event)[0]:
#             print("over %s" % curve.get_gid())

# import plotly.express as px
#
import plotly
from plotly import graph_objects as go
import scipy


fig = go.Figure()

fig.add_trace(go.Scatter(x=scan1, y=xx,
                    mode='lines+markers',
                    name='lines+markers'))

fig.write_html("original_TIC.html")

# fig, ax = plt.subplots(figsize=(6,1))
# ax.set_title("Original Chirp Signal: ")
# ax.plot(x= chirp_signal)


data = chirp_signal

waveletname1 = 'sym5'
waveletname= 'db28'

print(pywt.wavelist(kind='discrete'))

scales = np.arange(1, len(data))

fig, axarr = plt.subplots(nrows=5, ncols=2, figsize=(6,6))
for ii in range(5):
    (data, coeff_d) = pywt.dwt(data, wavelet= waveletname)

    axarr[ii, 0].plot(data, 'r')
    axarr[ii, 1].plot(coeff_d, 'g')
    axarr[ii, 0].set_ylabel("Level {}".format(ii + 1), fontsize=14, rotation=90)
    axarr[ii, 0].set_yticklabels([])
    if ii == 0:
        axarr[ii, 0].set_title("Approximation coefficients", fontsize=14)
        axarr[ii, 1].set_title("Detail coefficients", fontsize=14)
    axarr[ii, 1].set_yticklabels([])
plt.tight_layout()
plt.savefig('comparison_filter_step.svg', bbox_inches='tight', dpi= 1200)


def smoothing(data):
    for ii in range(3):
        (data, coeff_d) = pywt.dwt(data, 'db28')
        fig = go.Figure()

        fig.add_trace(go.Scatter(x=scan1, y=data,
                            mode='lines+markers',
                            name=f'smoothed spectra after_{ii}'))
        fig.write_html(f"smoothed_{ii}.html")
    return data

from scipy import signal
data= smoothing(xx)

peak_indexes = signal.argrelextrema(data, np.greater)
peak_indexes = peak_indexes[0]

# Find valleys(min).
valley_indexes = signal.argrelextrema(data, np.less)
valley_indexes = valley_indexes[0]

# Plot main graph.

fig = go.Figure()
fig.add_trace(go.Scatter(x=scan1, y=data,
                    mode='lines+markers',
                    name='Base_spectra'))


# Plot peaks.
peak_x = peak_indexes
peak_y = data[peak_indexes]
fig.add_trace(go.Scatter(x=peak_x, y=peak_y,
                    mode='lines+markers',
                    name='peaks_selected_top'))

# Plot valleys.
valley_x = valley_indexes
valley_y = data[valley_indexes]
fig.add_trace(go.Scatter(x=valley_x, y=valley_y,
                    mode='lines+markers',
                    name='peaks_selected_valley'))


# Save graph to file.
fig.write_html(f"picked_peaks.html")


df1= df.iloc[peak_indexes, :]

df1.to_csv('peaks.csv')


def peak_area(scan_array, intensity_array, start, stop):
    area = 0

    for i in range(start + 1, stop):
        x1 = scan_array[i - 1]
        y1 = intensity_array[i - 1]
        x2 = scan_array[i]
        y2 = intensity_array[i]
        area += (y1 * (x2 - x1)) + ((y2 - y1) * (x2 - x1) / 2.)

    return area


y= np.array(df['m/z'].tolist())
x= np.array(df['intensity'].tolist())
yy= np.array(df['rt'].tolist())
ss= np.array(df['scan'].tolist())


prominences, left_bases, right_bases = scipy.signal.peak_prominences(x, peak_indexes)

df2= pd.DataFrame()

df2['prominence']= prominences
df2['left_bases_scan']= left_bases

df2['right_bases_scan']= right_bases

df2['peaks_m/z']= y[peak_indexes]
df2['peaks_intensity']= x[peak_indexes]

df2.to_csv('properties.csv')



mass=[]

for i, z in zip(left_bases, right_bases):
    ass= []
    ass.append(y[i:z])
    mass.append(ass)

intensity= []


for i, t in zip(left_bases, right_bases):
    sity= []
    sity.append(x[i:t])
    intensity.append(sity)



time=[]

for i, z in zip(left_bases, right_bases):
    ass= []
    ass.append(yy[i:z])
    time.append(ass)



dfmass= pd.DataFrame()

dfmass['mass_range']= mass

dfintensity= pd.DataFrame()

dfintensity['intensity_range']= intensity

dfintensity['time_range']= time


df3= pd.concat([dfmass, dfintensity, df2], axis=1)

print(list(df3.columns))


df3['number_of_points']= df3['mass_range'].map(lambda x: len(x[0]))

width= peak_widths(x, np.array(peak_indexes), rel_height=1)

width= width[0]

df3['peak_width']= width

def row_area(row):
    cc = row['intensity_range']
    ck = row['intensity_range']
    cm = row['number_of_points']
    return peak_area(list(range(row['number_of_points'])), ck[0], 0, cm)



df3['area']= df3.apply(row_area, axis= 1)


df3.to_csv('peak_full_analysis.csv')
