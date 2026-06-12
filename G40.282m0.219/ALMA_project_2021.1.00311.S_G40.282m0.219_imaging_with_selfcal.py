import os

# the script was used with CASA version 6.6.1-17

# an example of how to run the script:
# casa-6.6.1-17-pipeline-2024.1.0.8/bin/mpicasa -n 14 casa-6.6.1-17-pipeline-2024.1.0.8/bin/casa --nogui --log2term --logfile "log.txt" -c "execfile('ALMA_project_2021.1.00311.S_G40.282m0.219_imaging_with_selfcal.py')" &> term_log.txt

cont_channels = '25:216.6550264819~216.6726046069GHz;216.6921358569~216.7927217944GHz;216.8542452319~216.9069796069GHz;216.9489717944~216.9704561694GHz;216.9763155444~216.9890108569GHz;216.9977999194~217.0114717944GHz;217.0212374194~217.0534639819GHz;217.1140108569~217.1267061694GHz;217.1374483569~217.1452608569GHz;217.1501436694~217.1569796069GHz;217.1686983569~217.1813936694GHz;217.1911592944~217.1999483569GHz;217.2175264819~217.2292452319GHz;217.2380342944~217.2634249194GHz;217.2761202319~217.2790499194GHz;217.2878389819~217.3649874194GHz;217.3757296069~217.3874483569GHz;217.4030733569~217.4157686694GHz;217.4255342944~217.4587374194GHz;217.4704561694~217.4870577319GHz;217.4987764819~217.5163546069GHz;217.5231905444~217.5446749194GHz;217.5700655444~217.5759249194GHz;217.5866671069~217.6032686694GHz;217.6198702319~217.6462374194GHz;217.6882296069~217.7077608569GHz;217.7194796069~217.7702608569GHz;217.8220186694~217.8298311694GHz;217.8327608569~217.8483858569GHz;217.8669405444~217.8903780444GHz;217.9304171069~217.9724092944GHz;218.0290499194~218.0446749194GHz;218.0847139819~218.0925264819GHz;218.1052217944~218.1091280444GHz;218.1198702319~218.1227999194GHz;218.1384249194~218.1442842944GHz;218.1628389819~218.1657686694GHz;218.2116671069~218.2253389819GHz;218.2331514819~218.2448702319GHz;218.2595186694~218.2634249194GHz;218.2770967944~218.2819796069GHz;218.3102999194~218.3376436694GHz;218.3806124194~218.3913546069GHz;218.4743624194~218.4860811694GHz;218.4958467944~218.5065889819GHz;218.5144014819~218.5212374194GHz,' + \
 '27:219.0294905713~219.0353499463GHz;219.0412093213~219.0470686963GHz;219.0509749463~219.0548811963GHz;219.0675765088~219.0802718213GHz;219.0832015088~219.1212874463GHz;219.1691390088~219.1857405713GHz;219.2521468213~219.2687483838GHz;219.3361311963~219.3488265088GHz;219.3898421338~219.4035140088GHz;219.4191390088~219.4318343213GHz;219.4513655713~219.4699202588GHz;219.4865218213~219.5050765088GHz;219.5578108838~219.5665999463GHz;219.5763655713~219.5880843213GHz;219.5988265088~219.6085921338GHz;219.7169905713~219.7287093213GHz;219.7394515088~219.7472640088GHz;219.7501936963~219.7609358838GHz;219.7999983838~219.8419905713GHz;219.8498030713~219.8595686963GHz;219.8664046338~219.8722640088GHz;219.9708968213~219.9943343213GHz;220.0177718213~220.0363265088GHz;220.0607405713~220.0900374463GHz;220.1915999463~220.2140608838GHz;220.2365218213~220.2482405713GHz;220.4201155713~220.4279280713GHz;220.4728499463~220.4982405713GHz;220.5783186963~220.5871077588GHz;220.7658186963~220.7716780713GHz;220.8751936963~220.8771468213GHz;220.8878890088~220.8898421338GHz,' + \
 '29:230.0873039225~230.1009757975GHz;230.1166007975~230.1253898600GHz;230.1556632975~230.1722648600GHz;230.1849601725~230.1947257975GHz;230.2113273600~230.2201164225GHz;230.2962882975~230.3070304850GHz;230.3138664225~230.3187492350GHz;230.3499992350~230.3626945475GHz;230.3714836100~230.3900382975GHz;230.5707023600~230.5785148600GHz;230.6029289225~230.6126945475GHz;230.6214836100~230.6292961100GHz;230.6546867350~230.6644523600GHz;230.6859367350~230.6947257975GHz;230.6976554850~230.7015617350GHz;230.7132804850~230.7386711100GHz;230.7464836100~230.7533195475GHz;230.7816398600~230.7933586100GHz;230.8050773600~230.8412101725GHz;230.8675773600~230.8949211100GHz;230.9115226725~230.9173820475GHz;230.9447257975~230.9544914225GHz;230.9681632975~230.9828117350GHz;230.9955070475~231.0091789225GHz;231.0511711100~231.0804679850GHz;231.1244132975~231.1322257975GHz;231.2269523600~231.2328117350GHz;231.2650382975~231.2718742350GHz;231.2884757975~231.2943351725GHz;231.3939445475~231.4046867350GHz;231.4085929850~231.4144523600GHz;231.4408195475~231.4496086100GHz;231.4593742350~231.4662101725GHz;231.4964836100~231.5160148600GHz;231.5404289225~231.5492179850GHz;231.5902336100~231.6205070475GHz;231.6605461100~231.6800773600GHz;231.8460929850~231.8548820475GHz;231.8773429850~231.8851554850GHz,' + \
 '31:232.1062011216~232.1101073716GHz;232.4782714341~232.5349120591GHz;232.5495604966~232.5573729966GHz;232.5700683091~232.5808104966GHz;232.6315917466~232.6374511216GHz;232.7116698716~232.7145995591GHz;232.7360839341~232.7399901841GHz;232.8581542466~232.8757323716GHz;233.2038573716~233.2116698716GHz;233.3386229966~233.3503417466GHz;233.3972167466~233.4060058091GHz'

def make_cont_image(ivis, imname, ithresh='0.00178Jy'):
    os.system('rm -rf '+imname+'*')
    tclean(vis=ivis, field='G40.282-0.219', spw=cont_channels, intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='corrected', imagename=imname, imsize=[256, 256], cell=['0.092arcsec'], phasecenter='ICRS 19:05:41.2200 +006.26.12.700', stokes='I', specmode='mfs', nchan=-1, outframe='LSRK', perchanweightdensity=False, gridder='standard', mosweight=False, usepointing=False, pblimit=0.2, deconvolver='asp', restoration=False, restoringbeam='common', pbcor=False, weighting='briggs', robust=0.5, npixels=0, niter=0, threshold='0.0mJy', nsigma=0.0, interactive=False, mask="circle[[286.421733434deg, 6.437012912deg], 7.1arcsec]", sidelobethreshold=2.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, minbeamfrac=0.3, growiterations=75, dogrowprune=True, minpercentchange=1.0, fastnoise=False, savemodel='modelcolumn', largestscale=100, parallel=True)
    tclean(vis=ivis, field='G40.282-0.219', spw=cont_channels, intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='corrected', imagename=imname, imsize=[256, 256], cell=['0.092arcsec'], phasecenter='ICRS 19:05:41.2200 +006.26.12.700', stokes='I', specmode='mfs', nchan=-1, outframe='LSRK', perchanweightdensity=False, gridder='standard', mosweight=False, usepointing=False, pblimit=0.2, deconvolver='asp', restoration=True, restoringbeam='common', pbcor=True, weighting='briggs', robust=0.5, npixels=0, niter=7000000, threshold=ithresh, nsigma=0.0, interactive=False, mask="circle[[286.421733434deg, 6.437012912deg], 7.1arcsec]", sidelobethreshold=2.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, minbeamfrac=0.3, growiterations=75, dogrowprune=True, minpercentchange=1.0, fastnoise=False, restart=True, savemodel='modelcolumn', calcres=False, calcpsf=False, largestscale=100, parallel=True)
    exportfits(imname+".image.pbcor", fitsimage=imname+".image.pbcor.fits", overwrite=True)

def make_cubes(ivis, imname, ithresh='0.00178Jy', ispw="25"):
    os.system('rm -rf '+imname+'*')
    inchan = -1
    istart = ''
    iwidth = ''
    if ispw == "25":
        inchan = 1918
        istart = '216.6306922414GHz'
        iwidth = '0.9764842MHz'
    if ispw == "27":
        inchan = 1917
        istart = '219.0000020321GHz'
        iwidth = '0.9764845MHz'
    if ispw == "29":
        inchan = 1918
        istart = '230.0467839682GHz'
        iwidth = '0.9764862MHz'
    if ispw == "31":
        inchan = 1918
        istart = '231.8662539342GHz'
        iwidth = '0.9764865MHz'
    tclean(vis=ivis, field='G40.282-0.219', spw=ispw, intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='corrected', imagename=imname, imsize=[392, 392], cell=['0.092arcsec'], phasecenter='ICRS 19:05:41.2200 +006.26.12.700', stokes='I', specmode='cube', nchan=inchan, start=istart, width=iwidth, outframe='LSRK', perchanweightdensity=True, gridder='standard', mosweight=False, usepointing=False, pblimit=0.2, deconvolver='asp', restoration=False, restoringbeam='common', pbcor=False, weighting='briggsbwtaper', robust=0.5, npixels=0, niter=0, threshold='0.0mJy', nsigma=0.0, interactive=False, mask="circle[[286.421673582deg, 6.437055333deg], 11.7arcsec]", sidelobethreshold=2.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=7.0, minbeamfrac=0.3, growiterations=50, dogrowprune=True, minpercentchange=1.0, fastnoise=False, savemodel='modelcolumn', largestscale=100, parallel=True)
    tclean(vis=ivis, field='G40.282-0.219', spw=ispw, intent='OBSERVE_TARGET#ON_SOURCE', datacolumn='corrected', imagename=imname, imsize=[392, 392], cell=['0.092arcsec'], phasecenter='ICRS 19:05:41.2200 +006.26.12.700', stokes='I', specmode='cube', nchan=inchan, start=istart, width=iwidth, outframe='LSRK', perchanweightdensity=True, gridder='standard', mosweight=False, usepointing=False, pblimit=0.2, deconvolver='asp', restoration=True, restoringbeam='common', pbcor=True, weighting='briggsbwtaper', robust=0.5, npixels=0, niter=4000000, threshold=ithresh, nsigma=0.0, interactive=False, mask="circle[[286.421673582deg, 6.437055333deg], 11.7arcsec]", sidelobethreshold=2.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=7.0, minbeamfrac=0.3, growiterations=50, dogrowprune=True, minpercentchange=1.0, fastnoise=False, restart=True, savemodel='modelcolumn', calcres=False, calcpsf=False, largestscale=100, parallel=True)
    exportfits(imname+".image.pbcor", fitsimage=imname+".image.pbcor.fits", overwrite=True)


visfile = 'uid___A002_Xf9cb82_X2437.ms.split.cal'

'''
h_init()
try:
    hifa_importdata(vis=visfile, dbservice=False, pipelinemode="automatic")
    hifa_flagtargets(pipelinemode="automatic")
    hif_uvcontsub(pipelinemode="automatic")
finally:
    h_save()
'''

imname ='G40.282-0.219_cont_self_cal_0'
make_cont_image(visfile, imname, ithresh='0.00178Jy')

os.system('rm -rf '+ visfile + '.phase_int.tb1')
gaincal(vis=visfile,
        caltable = visfile + '.phase_int.tb1',
        field = 'G40.282-0.219',
        spw=cont_channels,
        solint='inf',
        calmode='p',
        refant='DA63, DV08, DA52',
        gaintype='G')
applycal(vis=visfile,
         field = 'G40.282-0.219',
         gaintable=[visfile + '.phase_int.tb1'],
         flagbackup = False)
os.system('rm -rf '+ visfile + "1*")
split(vis=visfile,
      outputvis=visfile + "1",
      datacolumn='corrected')

imname ='G40.282-0.219_cont_self_cal_1'
make_cont_image(visfile + "1", imname, ithresh='0.0012Jy')

os.system('rm -rf '+ visfile + '.phase_int.tb2')
gaincal(vis=visfile + "1",
        caltable = visfile + '.phase_int.tb2',
        field = 'G40.282-0.219',
        spw=cont_channels,
        solint='int',
        calmode='p',
        refant='DA63, DV08, DA52',
        gaintype='G')
applycal(vis=visfile + "1",
         field = 'G40.282-0.219',
         gaintable=[visfile + '.phase_int.tb2'],
         flagbackup = False)
os.system('rm -rf '+ visfile + "2*")
split(vis=visfile + "1",
      outputvis=visfile + "2",
      datacolumn='corrected')

imname ='G40.282-0.219_cont_self_cal_2'
make_cont_image(visfile + "2", imname, ithresh='0.0010Jy')

visfileline = 'uid___A002_Xf9cb82_X2437_line.ms'

applycal(vis=visfileline,
         field = 'G40.282-0.219',
         gaintable=[visfile + '.phase_int.tb1'],
         flagbackup = False)
split(vis=visfileline,
      outputvis=visfileline + "1",
      datacolumn='corrected')

applycal(vis=visfileline + "1",
         field = 'G40.282-0.219',
         gaintable=[visfile + '.phase_int.tb2'],
         flagbackup = False)
split(vis=visfileline + "1",
      outputvis=visfileline + "2",
      datacolumn='corrected')

visfileline = visfileline + "2"

imname ='G40.282-0.219_cube_self_cal'
make_cubes(visfileline, imname + '_spw25', ithresh='0.0199Jy', ispw="25")
make_cubes(visfileline, imname + '_spw27', ithresh='0.0239Jy', ispw="27")
make_cubes(visfileline, imname + '_spw29', ithresh='0.023Jy', ispw="29")
make_cubes(visfileline, imname + '_spw31', ithresh='0.0195Jy', ispw="31")
