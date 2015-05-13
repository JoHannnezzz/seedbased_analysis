# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 11:49:45 2015

@author: jgolchert
"""
#Run in terminal before staring script
#export AFNI_1D_ZERO_TEXT=YES

import os                                    # system functions\n",
import nipype.interfaces.freesurfer as fs    # freesurfer\n",
import nipype.interfaces.io as nio           # i/o routines\n",
import nipype.interfaces.utility as util     # utility\n",
import nipype.pipeline.engine as pe          # pypeline engine\n",
import nipype.interfaces.fsl as fsl
import nipype.interfaces.afni as afni

#specify ROI-coordinates
rois = []  

#test Rois
#rois.append((30,-21,-18))
rois.append((2,-52,28))
#rois.append((30,-21,-18))


#DMN_ROIs (Andrews-Hanna, 2010; Neuron: Regions are 8mm spheres, we use 4mm radius for now

#PCC - aMPFCcore
#rois.append((-6,52,-2)) #aMPFC
#rois.append((6,52,-2))
#rois.append((-8,-56,26)) #PCC
#rois.append((8,-56,26))

#dMPFC subsytsem
#rois.append((0,52,26)) #dMPFC
#rois.append((-54,-54,28)) TPJ
#rois.append((54,-54,28))
#rois.append((-60,-24,-18)) LTC
#rois.append((60,-24,-18))
#rois.append((-50,14,-40)) Temporal Pole
#rois.append((50,14,-40))

#MTL subsystem
#rois.append((0,26,-18)) #vMPFC
#rois.append((-44,-74,32)) #pIPL
#rois.append((44,-74,32))
#rois.append((-14,-52,8)) retrosplenial cortex
#rois.append((14,-52,8))
#rois.append((-28,-40,-12)) parahipp. cortex
#rois.append((28,-40,-12))
#rois.append((-22,-20,-26)) Hippocamp.
#rois.append((22,-20,-26)

#subjects=['00796']
# read in subjects
sublist = '/scr/ilz2/LEMON_LSD/all_lsd_rest1a.txt'
with open(sublist, 'r') as f:
    subjects = [line.strip() for line in f]
print subjects

sess = ['rest1a','rest1b','rest2a','rest2b']
workingdir = "/home/raid1/margulies/working_dir/"
resultsdir = "/home/raid1/margulies/results_dir/"

wf = pe.Workflow(name="main")
wf.base_dir = workingdir
wf.config['execution']['crashdump_dir'] = wf.base_dir + "crash_files"

roi_infosource = pe.Node(util.IdentityInterface(fields=['roi']), name="roi_infosource")
roi_infosource.iterables = ('roi', rois) 

subjects_infosource = pe.Node(interface=util.IdentityInterface(fields=['subject_id']), name="subjects_infosource")
subjects_infosource.iterables = ("subject_id", subjects)

sess_infosource = pe.Node(interface=util.IdentityInterface(fields=['sess_id']), name="sess_infosource")
sess_infosource.iterables = ("sess_id", sess)

ds = pe.Node(nio.DataSink(), name="datasink")
#ds.run_without_submitting = True",
ds.inputs.base_directory = resultsdir

datasource = pe.Node(nio.DataGrabber(infields=['subject_id', 'sess_id'], outfields = ['EPI_bandpassed']), name="datasource") #grabs data
datasource.inputs.base_directory = '/scr/ilz2/LEMON_LSD/'
datasource.inputs.template = '%s/preprocessed/lsd_resting/%s/rest_preprocessed2mni.nii.gz'
datasource.inputs.template_args['EPI_bandpassed'] = [['subject_id', 'sess_id']] 
datasource.inputs.sort_filelist = True
wf.connect(subjects_infosource, "subject_id", datasource, "subject_id")
wf.connect(sess_infosource, "sess_id", datasource, "sess_id")

automask = pe.Node(interface=afni.Automask(), name='automask')
automask.inputs.dilate = 1
automask.inputs.outputtype = "NIFTI_GZ"
wf.connect(datasource, 'EPI_bandpassed', automask, 'in_file')
wf.connect(automask, 'out_file', ds, '@automask')

#extract rois with spheres
sphere = pe.Node(afni.Calc(), name="sphere")
sphere.inputs.in_file_a = fsl.Info.standard_image('/usr/share/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz')
sphere.inputs.outputtype='NIFTI_GZ'

def roi2exp(coord): 
    radius = 4
    return "step((%d*%d)-(x+%d)*(x+%d)-(y+%d)*(y+%d)-(z+%d)*(z+%d))" %(radius, radius, coord[0], coord[0], coord[1], coord[1], -coord[2], -coord[2])

def roi2name(coord):
    return 'roi_sphere_%s_%s_%s.nii.gz'%(str(coord[0]), str(coord[1]), str(coord[2]))
    
wf.connect(roi_infosource, ("roi", roi2exp), sphere, "expr")   
wf.connect(roi_infosource, ("roi", roi2name), sphere,"out_file") 
wf.connect(sphere, "out_file", ds, "@sphere")

extract_timeseries = pe.Node(afni.Maskave(), name="extract_timeseries") 
extract_timeseries.inputs.quiet = True
wf.connect(sphere, "out_file", extract_timeseries, "mask") 
wf.connect(datasource, 'EPI_bandpassed', extract_timeseries, "in_file") 

correlation_map = pe.Node(afni.Fim(), name="correlation_map")
correlation_map.inputs.out = "Correlation"
correlation_map.inputs.outputtype = "NIFTI_GZ"
correlation_map.inputs.out_file = "corr_map.nii.gz"
correlation_map.inputs.fim_thr = 0.0001
wf.connect(extract_timeseries, "out_file", correlation_map, "ideal_file")
wf.connect(datasource, 'EPI_bandpassed', correlation_map, "in_file")
def add_arg_fim_mask(input):
    return ('-mask %s' % input)
wf.connect(automask, ('out_file', add_arg_fim_mask), correlation_map, 'args')

z_trans = pe.Node(interface=afni.Calc(), name='z_trans')
z_trans.inputs.expr = 'log((1+a)/(1-a))/2'
z_trans.inputs.outputtype = 'NIFTI_GZ'
wf.connect(correlation_map, "out_file", z_trans, "in_file_a")
wf.connect(z_trans, 'out_file', ds, "@z_trans")

# Smoothing now happens *after* correlation map is caluclated
smooth = pe.Node(fsl.maths.IsotropicSmooth(), name = "smooth")
smooth.inputs.fwhm = 6.0
smooth.inputs.out_file = "corr_map_smoothed.nii.gz"
def add_arg_smooth(input):
    return ('-mas %s' % input)
wf.connect(automask, ('out_file', add_arg_smooth), smooth, 'args')
wf.connect(z_trans, 'out_file', smooth, "in_file")
wf.connect(smooth, 'out_file', ds, "@smooth")

# take mean of correlation maps from 4 runs
def format_roi(roi_str):
    import string
    valid_chars = "-_.%s%s" % (string.ascii_letters, string.digits)
    return ''.join(c for c in str(roi_str).replace(",",".") if c in valid_chars)
    
func_datasource = pe.Node(nio.DataGrabber(infields=['subject_id', 'roi'], outfields = ['corr_maps']), name="func_datasource")
func_datasource.inputs.base_directory = workingdir
func_datasource.inputs.template = 'main/_roi_%s/_sess_id_*/_subject_id_%s/smooth/corr_map_smoothed.nii.gz'
func_datasource.inputs.template_args['corr_maps'] = [['roi', 'subject_id']]
func_datasource.inputs.sort_filelist = True
wf.connect(roi_infosource, ("roi", format_roi), func_datasource, "roi")  
wf.connect(subjects_infosource, 'subject_id', func_datasource, "subject_id")       

merge_conn_maps = pe.Node(fsl.Merge(dimension='t'), name="merge_conn_maps")
wf.connect(func_datasource, 'corr_maps', merge_conn_maps, 'in_files')

mean_conn_map = pe.Node(fsl.MeanImage(dimension="T"), name="mean_conn_map")
wf.connect(merge_conn_maps, 'merged_file', mean_conn_map, 'in_file')
wf.connect(mean_conn_map, 'out_file', ds, "mean_conn_maps")

std_conn_maps = pe.Node(fsl.ImageMaths(op_string = "-Tstd"),name="std_conn_maps")
wf.connect(merge_conn_maps, 'merged_file', std_conn_maps, 'in_file')
wf.connect(std_conn_maps, 'out_file', ds, "std_conn_maps")

#wf.write_graph(dotfilename='/scr/lahn2/jgolchert/LSD/MW_impulsivity/graph.dot', graph2use='colored', format ='pdf', simple_form=True)
wf.run()

