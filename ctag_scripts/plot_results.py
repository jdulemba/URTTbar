from argparse import ArgumentParser
import os, glob, sys
import URAnalysis.Utilities.prettyjson as prettyjson
from pdb import set_trace
import logging
import matplotlib.pyplot as plt
from styles import styles

parser = ArgumentParser()
args = parser.parse_args()

print 'CHECK TAGGERS TO MAKE SURE WHICH ONES ARE INCLUDED'

input_dir = '%s/plots/%s/ctageff/mass_discriminant/' % (os.environ['URA_PROJECT'],os.environ['jobid'])
out_dir = '%s/plots/%s/ctageff/' % (os.environ['URA_PROJECT'],os.environ['jobid'])

results = prettyjson.loads(open('%s/results.json' % input_dir).read())

#working_points = ['DeepCSVctagLoose', 'DeepCSVctagMedium', 'DeepCSVctagTight']
working_points = results.keys()
#set_trace()

    ## create lists for loose (L), medium (M), tight (T), and all wps to use in scatterplot
L_tag_means = []
M_tag_means = []
T_tag_means = []
L_tag_sf_upper = []
L_tag_sf_lower = []
M_tag_sf_upper = []
M_tag_sf_lower = []
T_tag_sf_upper = []
T_tag_sf_lower = []
L_tag_stat_upper = []
L_tag_stat_lower = []
M_tag_stat_upper = []
M_tag_stat_lower = []
T_tag_stat_upper = []
T_tag_stat_lower = []

tag_means = []
tag_sf_upper = []
tag_sf_lower = []
tag_stat_upper = []
tag_stat_lower = []
tagger_names = []

xaxis_title = []
for wp in sorted(working_points):
        # find tag name and wp (L, M, T)
    name = results[wp].keys()[0]
    if name[:-2] not in xaxis_title:
        xaxis_title.append(name[:-2])

        # assign vars from json
    mean = results[wp][name]['mean'] #scale factor mean
    if results[wp][name]['sf_upper'] <= 0:
        sf_upper = -1*results[wp][name]['sf_upper'] # total sf upper error
    else:
        sf_upper = results[wp][name]['sf_upper'] # total sf upper error
    if results[wp][name]['sf_lower'] <= 0:
        sf_lower = -1*results[wp][name]['sf_lower'] # total sf lower error changed to be positive
    else:
        sf_lower = results[wp][name]['sf_lower'] # total sf lower error changed to be positive
    if results[wp][name]['stat_upper'] <= 0:
        stat_upper = -1*results[wp][name]['stat_upper'] # upper error only from statistics
    else:
        stat_upper = results[wp][name]['stat_upper'] # upper error only from statistics
    if results[wp][name]['stat_lower'] <= 0:
        stat_lower = -1*results[wp][name]['stat_lower'] # lower error only from statistics
    else:
        stat_lower = results[wp][name]['stat_lower'] # lower error only from statistics

    #set_trace()
    if ' L' in name:
        L_tag_means.append(mean)
        L_tag_sf_upper.append(sf_upper)
        L_tag_sf_lower.append(sf_lower)
        L_tag_stat_upper.append(stat_upper)
        L_tag_stat_lower.append(stat_lower)
    elif ' M' in name:
        M_tag_means.append(mean)
        M_tag_sf_upper.append(sf_upper)
        M_tag_sf_lower.append(sf_lower)
        M_tag_stat_upper.append(stat_upper)
        M_tag_stat_lower.append(stat_lower)
    elif ' T' in name:
        T_tag_means.append(mean)
        T_tag_sf_upper.append(sf_upper)
        T_tag_sf_lower.append(sf_lower)
        T_tag_stat_upper.append(stat_upper)
        T_tag_stat_lower.append(stat_lower)
    else:
        logging.error('Not a valid working point.')
        sys.exit()
        
    tag_means.append(mean)
    tag_sf_upper.append(sf_upper)
    tag_sf_lower.append(sf_lower)
    tag_stat_upper.append(stat_upper)
    tag_stat_lower.append(stat_lower)
    tagger_names.append(name)

## create arrays containing sum of sf+stat uncertainties
tag_total_lower = [tag_stat_lower[i]+tag_sf_lower[i] for i in range(len(tag_stat_lower))]
tag_total_upper = [tag_stat_upper[i]+tag_sf_upper[i] for i in range(len(tag_stat_upper))]

plot_title = styles['data*']['name']

fig1 = plt.figure()
plt.title(plot_title)
plt.ylabel('Charm SF values')
plt.xlabel('Taggers')
#plt.axis([0,2,0.5,1.5])

#set_trace()
xaxis = list(range(len(L_tag_means)))
plt.errorbar(xaxis, L_tag_means, yerr=[L_tag_stat_lower, L_tag_stat_upper], color='r', label='L Stat.', fmt='ko', elinewidth=3, capthick=2)
plt.errorbar(xaxis, M_tag_means, yerr=[M_tag_stat_lower, M_tag_stat_upper], color='b', label='M Stat.', fmt='ko', elinewidth=3, capthick=2)
plt.errorbar(xaxis, T_tag_means, yerr=[T_tag_stat_lower, T_tag_stat_upper], color='k', label='T Stat.', fmt='ko', elinewidth=3, capthick=2)
plt.legend(numpoints=1)
plt.xlim(min(xaxis)-1, max(xaxis)+2)
plt.xticks(xaxis, xaxis_title, rotation=90)
plt.grid()
plt.tight_layout()
#set_trace()
fname1 = '%s%s_SF_stat_uncs_results' % (out_dir, plot_title)
fig1.savefig('%s.png' % fname1)

fig2 = plt.figure()
plt.title(plot_title)
plt.ylabel('Charm SF values')
plt.xlabel('Taggers')
#plt.axis([0,2,0.5,1.5])

#set_trace()
xaxis = list(range(len(L_tag_means)))
plt.errorbar(xaxis, L_tag_means, yerr=[L_tag_sf_lower, L_tag_sf_upper], color='r', label='L Stat+Sys', fmt='o', elinewidth=3, capthick=2)
plt.errorbar(xaxis, M_tag_means, yerr=[M_tag_sf_lower, M_tag_sf_upper], color='b', label='M Stat+Sys', fmt='o', elinewidth=3, capthick=2)
plt.errorbar(xaxis, T_tag_means, yerr=[T_tag_sf_lower, T_tag_sf_upper], color='k', label='T Stat+Sys', fmt='o', elinewidth=3, capthick=2)
plt.legend(numpoints=1)
plt.xlim(min(xaxis)-1, max(xaxis)+2)
plt.xticks(xaxis, xaxis_title, rotation=90)
plt.grid()
plt.tight_layout()
#set_trace()
fname2 = '%s%s_SF_sys_uncs_results' % (out_dir, plot_title)
fig2.savefig('%s.png' % fname2)

#set_trace()
fig3 = plt.figure()
plt.title(plot_title)
plt.ylabel('Charm SF values')
plt.xlabel('Taggers')
#plt.axis([0,2,0.5,1.5])
xaxis = list(range(len(tag_means)))

plt.errorbar(xaxis, tag_means, yerr=[tag_stat_lower, tag_stat_upper], color='k', fmt='o', label='Stat.', elinewidth=3, capthick=2)
plt.legend(numpoints=1)
plt.xlim(min(xaxis)-1, max(xaxis)+2)
plt.xticks(xaxis, tagger_names, rotation=90)
plt.grid()
plt.tight_layout()
fname3 = '%s%s_SF_stat_uncs_indiv_results' % (out_dir, plot_title)
fig3.savefig('%s.png' % fname3)

fig4 = plt.figure()
plt.title(plot_title)
plt.ylabel('Charm SF values')
plt.xlabel('Taggers')
#plt.axis([0,2,0.5,1.5])
xaxis = list(range(len(tag_means)))

#plt.errorbar(xaxis, tag_means, yerr=[tag_sf_lower, tag_sf_upper], color='k', fmt='o', label='Stat+Sys', ecolor='b', elinewidth=3, capsize=10, capthick=5)
plt.errorbar(xaxis, tag_means, yerr=[tag_sf_lower, tag_sf_upper], color='k', fmt='ko', label='Stat+Sys', ecolor='r', elinewidth=10, alpha=0.5)
plt.legend(numpoints=1)
plt.xlim(min(xaxis)-1, max(xaxis)+2)
plt.xticks(xaxis, tagger_names, rotation=90)
plt.grid()
plt.tight_layout()
fname4 = '%s%s_SF_sys_uncs_indiv_results' % (out_dir, plot_title)
fig4.savefig('%s.png' % fname4)

fig5 = plt.figure()
plt.title(plot_title)
plt.ylabel('Charm SF values')
plt.xlabel('Taggers')
#plt.axis([0,2,0.5,1.5])
xaxis = list(range(len(tag_means)))

#plt.errorbar(xaxis, tag_means, yerr=[tag_total_lower, tag_total_upper], color='k', fmt='ko', label='Sys+Stat.', ecolor='b', elinewidth=10)
plt.errorbar(xaxis, tag_means, yerr=[tag_stat_lower, tag_stat_upper], color='k', fmt='ko', label='Stat.', ecolor='k', capsize=4)
plt.errorbar(xaxis, tag_means, yerr=[tag_sf_lower, tag_sf_upper], color='k', fmt='ko', label='Stat+Sys', ecolor='r', elinewidth=10, alpha=0.5)
plt.legend(numpoints=1)
plt.xlim(min(xaxis)-1, max(xaxis)+2)
plt.xticks(xaxis, tagger_names, rotation=90)
plt.ylim(0, 1.2)
plt.grid()
plt.tight_layout()
fname5 = '%s%s_SF_both_uncs_indiv_results' % (out_dir, plot_title)
fig5.savefig('%s.png' % fname5)

#fig6 = plt.figure()
#plt.title(plot_title)
#plt.ylabel('Charm SF values')
#plt.xlabel('Taggers')
#dEB.method2(xaxis, tag_means, [tag_stat_lower, tag_stat_upper], [tag_sf_lower, tag_sf_upper])
##plt.legend(numpoints=1)
#plt.xlim(min(xaxis)-1, max(xaxis)+2)
#plt.xticks(xaxis, tagger_names, rotation=90)
#plt.grid()
#plt.tight_layout()
#fname6 = '%s%s_SF_both_uncs_indiv_results' % (out_dir, plot_title)
#fig6.savefig('%s.png' % fname6)



#print '\n%s.png\n%s.png\n%s.png\n  have been created.' % (fname1, fname2, fname5)
print '\n%s.png\n%s.png\n%s.png\n%s.png\n%s.png\n  have been created.' % (fname1, fname2, fname3, fname4, fname5)
