#from xvfbwrapper import Xvfb

#vdisplay = Xvfb()
#vdisplay.start()

# the script was used with CASA version 6.6.1-17

h_init()
try:
    hifa_importdata(vis=['uid___A002_Xfcf230_Xb54.ms', 'uid___A002_Xfcf230_X1583.ms'], session=['session_1', 'session_2'],dbservice=False,pipelinemode="automatic")
    hifa_imageprecheck(pipelinemode="automatic", parallel=True)
    hif_mstransform(pipelinemode="automatic")
    hifa_flagtargets(pipelinemode="automatic")
    hif_makeimlist(field='G25.650+1.049', specmode='mfs')
    hif_findcont(pipelinemode="automatic", parallel=True)
    hif_uvcontsub(pipelinemode="automatic")
    hif_makeimages(pipelinemode="automatic", parallel=True)
    hif_makeimlist(field='G25.650+1.049', specmode='cont')
    hif_makeimages(pipelinemode="automatic", parallel=True)
    hif_makeimlist(field='G25.650+1.049', hm_imsize=[540, 540], specmode='cube')
    hif_makeimages(pipelinemode="automatic", parallel=True)
    hif_selfcal(field='G25.650+1.049')
    hif_makeimlist(field='G25.650+1.049', specmode='cont', datatype='selfcal')
    hif_makeimages(parallel=True)
    hif_makeimlist(field='G25.650+1.049', hm_imsize=[540, 540], specmode='cube', datatype='selfcal')
    hif_makeimages(parallel=True)
    hifa_exportdata(imaging_products_only=True)
finally:
    h_save()
    #vdisplay.stop()
