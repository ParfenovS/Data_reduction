import os

# the script was used with CASA version 6.6.1-17

def make_clean(inpvis, imname, mfs=True, inpthreshold='0.0205Jy', inpcont_channels=["27"], inpdatacolumn='corrected'):
        os.system('rm -rf '+imname+'*')
        if mfs:
                tclean(vis=inpvis, field='"G12.89+0.49"', spw=inpcont_channels, intent='OBSERVE_TARGET#ON_SOURCE', datacolumn=inpdatacolumn, imagename=imname, imsize=[1440, 1440], cell=['0.035arcsec'], phasecenter='ICRS 18:11:51.2700 -017.31.28.500', stokes='I', specmode='mfs', perchanweightdensity=True, gridder='standard', mosweight=False, usepointing=False, pblimit=0.2, deconvolver='asp', mask=["circle[[272.96382592deg, -17.52478752deg], 5.0098arcsec]","circle[[272.96013352deg, -17.52456916deg], 2.1645arcsec]"], restoration=False, restoringbeam='common', pbcor=False, weighting='briggs', robust=0.5, npixels=0, niter=0, threshold='0.0mJy', nsigma=0.0, interactive=False, sidelobethreshold=2.5, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=7.0, minbeamfrac=0.3, growiterations=50, dogrowprune=True, minpercentchange=1.0, fastnoise=True, savemodel='modelcolumn', parallel=False)
                tclean(vis=inpvis, field='"G12.89+0.49"', spw=inpcont_channels, intent='OBSERVE_TARGET#ON_SOURCE', datacolumn=inpdatacolumn, imagename=imname, imsize=[1440, 1440], cell=['0.035arcsec'], phasecenter='ICRS 18:11:51.2700 -017.31.28.500', stokes='I', specmode='mfs', perchanweightdensity=True, gridder='standard', mosweight=False, usepointing=False, pblimit=0.2, deconvolver='asp', mask=["circle[[272.96382592deg, -17.52478752deg], 5.0098arcsec]","circle[[272.96013352deg, -17.52456916deg], 2.1645arcsec]"], restoration=True, restoringbeam='common', pbcor=True, weighting='briggs', robust=0.5, npixels=0, niter=9000000, threshold=inpthreshold, nsigma=0.0, interactive=False, sidelobethreshold=2.5, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=7.0, minbeamfrac=0.3, growiterations=50, dogrowprune=True, minpercentchange=1.0, fastnoise=True, restart=True, savemodel='modelcolumn', calcres=False, calcpsf=False, parallel=False)
        else:
                if inpcont_channels == ["25"]:
                        tclean(vis=inpvis, field='"G12.89+0.49"', spw=inpcont_channels, intent='OBSERVE_TARGET#ON_SOURCE', datacolumn=inpdatacolumn, imagename=imname, imsize=[540, 540], cell=['0.051arcsec'], phasecenter='ICRS 18:11:51.2700 -017.31.28.500', stokes='I', specmode='cube', nchan=1917, start='216.6685544217GHz', width='0.9765997MHz', outframe='LSRK', perchanweightdensity=True, gridder='standard', mosweight=False, usepointing=False, pblimit=0.2, deconvolver='asp', mask="circle[[272.96453755deg, -17.52379629deg], 9.0938arcsec]", restoration=False, restoringbeam='common', pbcor=False, weighting='briggsbwtaper', robust=0.5, npixels=0, niter=0, threshold='0.0mJy', nsigma=0.0, interactive=False, sidelobethreshold=2.5, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=7.0, minbeamfrac=0.3, growiterations=50, dogrowprune=True, minpercentchange=1.0, fastnoise=True, savemodel='modelcolumn', parallel=False)
                        tclean(vis=inpvis, field='"G12.89+0.49"', spw=inpcont_channels, intent='OBSERVE_TARGET#ON_SOURCE', datacolumn=inpdatacolumn, imagename=imname, imsize=[540, 540], cell=['0.051arcsec'], phasecenter='ICRS 18:11:51.2700 -017.31.28.500', stokes='I', specmode='cube', nchan=1917, start='216.6685544217GHz', width='0.9765997MHz', outframe='LSRK', perchanweightdensity=True, gridder='standard', mosweight=False, usepointing=False, pblimit=0.2, deconvolver='asp', mask="circle[[272.96453755deg, -17.52379629deg], 9.0938arcsec]", restoration=True, restoringbeam='common', pbcor=True, weighting='briggsbwtaper', robust=0.5, npixels=0, niter=9000000, threshold=inpthreshold, nsigma=0.0, interactive=False, sidelobethreshold=2.5, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=7.0, minbeamfrac=0.3, growiterations=50, dogrowprune=True, minpercentchange=1.0, fastnoise=True, restart=True, savemodel='modelcolumn', calcres=False, calcpsf=False, parallel=False)
                elif inpcont_channels == ["27"]:
                        tclean(vis=inpvis, field='"G12.89+0.49"', spw=inpcont_channels, intent='OBSERVE_TARGET#ON_SOURCE', datacolumn=inpdatacolumn, imagename=imname, imsize=[540, 540], cell=['0.051arcsec'], phasecenter='ICRS 18:11:51.2700 -017.31.28.500', stokes='I', specmode='cube', nchan=1918, start='219.0382779633GHz', width='0.9765998MHz', outframe='LSRK', perchanweightdensity=True, gridder='standard', mosweight=False, usepointing=False, pblimit=0.2, deconvolver='asp', mask="circle[[272.96453755deg, -17.52379629deg], 9.0938arcsec]", restoration=False, restoringbeam='common', pbcor=False, weighting='briggsbwtaper', robust=0.5, npixels=0, niter=0, threshold='0.0mJy', nsigma=0.0, interactive=False, sidelobethreshold=2.5, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=7.0, minbeamfrac=0.3, growiterations=50, dogrowprune=True, minpercentchange=1.0, fastnoise=True, savemodel='modelcolumn', parallel=False)
                        tclean(vis=inpvis, field='"G12.89+0.49"', spw=inpcont_channels, intent='OBSERVE_TARGET#ON_SOURCE', datacolumn=inpdatacolumn, imagename=imname, imsize=[540, 540], cell=['0.051arcsec'], phasecenter='ICRS 18:11:51.2700 -017.31.28.500', stokes='I', specmode='cube', nchan=1918, start='219.0382779633GHz', width='0.9765998MHz', outframe='LSRK', perchanweightdensity=True, gridder='standard', mosweight=False, usepointing=False, pblimit=0.2, deconvolver='asp', mask="circle[[272.96453755deg, -17.52379629deg], 9.0938arcsec]", restoration=True, restoringbeam='common', pbcor=True, weighting='briggsbwtaper', robust=0.5, npixels=0, niter=9000000, threshold=inpthreshold, nsigma=0.0, interactive=False, sidelobethreshold=2.5, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=7.0, minbeamfrac=0.3, growiterations=50, dogrowprune=True, minpercentchange=1.0, fastnoise=True, restart=True, savemodel='modelcolumn', calcres=False, calcpsf=False, parallel=False)
                elif inpcont_channels == ["29"]:
                        tclean(vis=inpvis, field='"G12.89+0.49"', spw=inpcont_channels, intent='OBSERVE_TARGET#ON_SOURCE', datacolumn=inpdatacolumn, imagename=imname, imsize=[540, 540], cell=['0.051arcsec'], phasecenter='ICRS 18:11:51.2700 -017.31.28.500', stokes='I', specmode='cube', nchan=1917, start='230.0869889909GHz', width='0.9766004MHz', outframe='LSRK', perchanweightdensity=True, gridder='standard', mosweight=False, usepointing=False, pblimit=0.2, deconvolver='asp', mask="circle[[272.96453755deg, -17.52379629deg], 9.0938arcsec]", restoration=False, restoringbeam='common', pbcor=False, weighting='briggsbwtaper', robust=0.5, npixels=0, niter=0, threshold='0.0mJy', nsigma=0.0, interactive=False, sidelobethreshold=2.5, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=7.0, minbeamfrac=0.3, growiterations=50, dogrowprune=True, minpercentchange=1.0, fastnoise=True, savemodel='modelcolumn', parallel=False)
                        tclean(vis=inpvis, field='"G12.89+0.49"', spw=inpcont_channels, intent='OBSERVE_TARGET#ON_SOURCE', datacolumn=inpdatacolumn, imagename=imname, imsize=[540, 540], cell=['0.051arcsec'], phasecenter='ICRS 18:11:51.2700 -017.31.28.500', stokes='I', specmode='cube', nchan=1917, start='230.0869889909GHz', width='0.9766004MHz', outframe='LSRK', perchanweightdensity=True, gridder='standard', mosweight=False, usepointing=False, pblimit=0.2, deconvolver='asp', mask="circle[[272.96453755deg, -17.52379629deg], 9.0938arcsec]", restoration=True, restoringbeam='common', pbcor=True, weighting='briggsbwtaper', robust=0.5, npixels=0, niter=9000000, threshold=inpthreshold, nsigma=0.0, interactive=False, sidelobethreshold=2.5, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=7.0, minbeamfrac=0.3, growiterations=50, dogrowprune=True, minpercentchange=1.0, fastnoise=True, restart=True, savemodel='modelcolumn', calcres=False, calcpsf=False, parallel=False)
                elif inpcont_channels == ["31"]:
                        tclean(vis=inpvis, field='"G12.89+0.49"', spw=inpcont_channels, intent='OBSERVE_TARGET#ON_SOURCE', datacolumn=inpdatacolumn, imagename=imname, imsize=[540, 540], cell=['0.051arcsec'], phasecenter='ICRS 18:11:51.2700 -017.31.28.500', stokes='I', specmode='cube', nchan=1917, start='231.9067766895GHz', width='0.9766005MHz', outframe='LSRK', perchanweightdensity=True, gridder='standard', mosweight=False, usepointing=False, pblimit=0.2, deconvolver='asp', mask="circle[[272.96453755deg, -17.52379629deg], 9.0938arcsec]", restoration=False, restoringbeam='common', pbcor=False, weighting='briggsbwtaper', robust=0.5, npixels=0, niter=0, threshold='0.0mJy', nsigma=0.0, interactive=False, sidelobethreshold=2.5, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=7.0, minbeamfrac=0.3, growiterations=50, dogrowprune=True, minpercentchange=1.0, fastnoise=True, savemodel='modelcolumn', parallel=False)
                        tclean(vis=inpvis, field='"G12.89+0.49"', spw=inpcont_channels, intent='OBSERVE_TARGET#ON_SOURCE', datacolumn=inpdatacolumn, imagename=imname, imsize=[540, 540], cell=['0.051arcsec'], phasecenter='ICRS 18:11:51.2700 -017.31.28.500', stokes='I', specmode='cube', nchan=1917, start='231.9067766895GHz', width='0.9766005MHz', outframe='LSRK', perchanweightdensity=True, gridder='standard', mosweight=False, usepointing=False, pblimit=0.2, deconvolver='asp', mask="circle[[272.96453755deg, -17.52379629deg], 9.0938arcsec]", restoration=True, restoringbeam='common', pbcor=True, weighting='briggsbwtaper', robust=0.5, npixels=0, niter=9000000, threshold=inpthreshold, nsigma=0.0, interactive=False, sidelobethreshold=2.5, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=7.0, minbeamfrac=0.3, growiterations=50, dogrowprune=True, minpercentchange=1.0, fastnoise=True, restart=True, savemodel='modelcolumn', calcres=False, calcpsf=False, parallel=False)



asdmname = 'uid___A002_Xfcc484_X8570'

msdata = asdmname + '.ms.split.cal'

####################################
# Continuum subtraction in UV plane
####################################

#cont_channels25 = "25:30~45;54~58;79~86;96~102;107~114;130~134;172~186;202~218;238~241;297~308;323~329;344~353;475~476;500~502;518~523;642~643;657~658;670~684;732~733;788~796;843~845;857~865;882~883;902~907;950~955;991~998;1002~1006;1051~1063;1120~1126;1188~1200;1231~1241;1266~1268;1301~1306;1316~1329;1396~1397;1437~1441;1454~1457;1491~1494;1761~1762;1779~1780;1854~1857;1900~1903"
cont_channels25 = "25:216.6923036256~216.706952246GHz;216.7157414182~216.719647717GHz;216.7401557856~216.7469918084GHz;216.7567575553~216.7626170035GHz;216.767499877~216.7743358998GHz;216.7899610949~216.7938673937GHz;216.830977232~216.8446492777GHz;216.8602744728~216.8758996679GHz;216.8954311618~216.8983608858GHz;216.9530490686~216.9637913903GHz;216.9784400107~216.9842994588GHz;216.9989480792~217.0077372515GHz;217.126879364~217.1278559387GHz;217.1512937313~217.1532468807GHz;217.1688720758~217.1737549493GHz;217.2899673378~217.2909439125GHz;217.3046159582~217.3055925328GHz;217.3173114292~217.3309834749GHz;217.3778590601~217.3788356348GHz;217.4325472429~217.4403598405GHz;217.4862588511~217.4882120005GHz;217.4999308968~217.5077434943GHz;217.5243452641~217.5253218388GHz;217.543876758~217.5487596314GHz;217.5907523432~217.5956352167GHz;217.6307919056~217.6376279285GHz;217.6415342273~217.645440526GHz;217.6893863872~217.7011052835GHz;217.756770041~217.7626294892GHz;217.8231771202~217.8348960165GHz;217.865169832~217.8749355789GHz;217.8993499462~217.9013030956GHz;217.9335300605~217.938412934GHz;217.9481786809~217.9608741519GHz;218.0263046563~218.027281231GHz;218.0663442187~218.0702505175GHz;218.0829459885~218.0858757126GHz;218.1190792522~218.1220089763GHz;218.3827544193~218.383730994GHz;218.4003327638~218.4013093385GHz;218.4735758658~218.4765055898GHz;218.5184983016~218.5214280257GHz"
#cont_channels27 = "27:50~57;98~99;189~190;231~236;316~326;399~405;524~526;641~643;710~714;769~773;788~796;803~811;842~843;883~886;925~936;961~981;1010~1012;1047~1067;1175~1191;1209~1212;1231~1233;1269~1271;1343~1350;1409~1415;1423~1431;1493~1497;1626~1629;1666~1667;1695~1697;1735~1736;1756~1766;1776~1780;1813~1822;1853~1861;1881~1890;1899~1910"
cont_channels27 = "27:219.0814979659~219.0883339894GHz;219.128373556~219.1293501308GHz;219.2172418621~219.2182184369GHz;219.2582580034~219.2631408774GHz;219.3412668608~219.3510326088GHz;219.4223225686~219.4281820174GHz;219.5443944178~219.5463475674GHz;219.6586536686~219.6606068181GHz;219.7260373293~219.7299436284GHz;219.7836552421~219.7875615412GHz;219.8022101631~219.8100227615GHz;219.816858785~219.8246713834GHz;219.8549452019~219.8559217767GHz;219.8949847685~219.8979144928GHz;219.9360009098~219.9467432325GHz;219.9711576023~219.9906890982GHz;220.0190097672~220.0209629167GHz;220.0551430345~220.0746745304GHz;220.180144608~220.1957698047GHz;220.213348151~220.2162778753GHz;220.2348327964~220.236785946GHz;220.2719426385~220.2738957881GHz;220.3442091732~220.3510451968GHz;220.4086631096~220.4145225583GHz;220.4223351567~220.430147755GHz;220.4906953922~220.4946016913GHz;220.6205798396~220.623509564GHz;220.6596428314~220.6606194062GHz;220.6879635004~220.6899166499GHz;220.7270264921~220.7280030669GHz;220.7475345627~220.7573003107GHz;220.7670660586~220.7709723578GHz;220.8031993259~220.8119884991GHz;220.8422623176~220.850074916GHz;220.8696064118~220.878395585GHz;220.8871847581~220.8979270808GHz"
#cont_channels29 = "29:4~10;15~23;37~40;78~87;133~142;216~232;274~283;302~305;343~345;354~361;466~472;580~595;652~661;684~690;727~730;761~777;820~826;833~835;849~852;881~894;907~910;959~963;993~998;1023~1028;1118~1122;1172~1179;1236~1241;1300~1305;1311~1316;1396~1400;1413~1416;1544~1545;1617~1620"
cont_channels29 = "29:230.08500543~230.0908648824GHz;230.0957477594~230.1035603626GHz;230.1172324182~230.1201621444GHz;230.1572720096~230.1660611882GHz;230.2109836567~230.2197728353GHz;230.2920394149~230.3076646213GHz;230.3486807882~230.3574699668GHz;230.3760248994~230.3789546256GHz;230.4160644908~230.4180176416GHz;230.4268068202~230.433642848GHz;230.5361832651~230.5420427175GHz;230.6475128607~230.6621614917GHz;230.7178262896~230.7266154682GHz;230.7490767024~230.7549361548GHz;230.7910694446~230.7939991708GHz;230.8242730082~230.8398982147GHz;230.8818909569~230.8877504093GHz;230.8945864371~230.8965395879GHz;230.9102116435~230.9131413697GHz;230.9414620563~230.9541575365GHz;230.9668530167~230.9697827429GHz;231.0176349376~231.0215412392GHz;231.0508385012~231.0557213782GHz;231.0801357632~231.0850186402GHz;231.1729104263~231.1768167279GHz;231.2256454979~231.2324815257GHz;231.2881463235~231.2930292005GHz;231.3506471492~231.3555300262GHz;231.3613894786~231.3662723556GHz;231.4443983876~231.4483046892GHz;231.4610001694~231.4639298956GHz;231.5889315469~231.5899081223GHz;231.6602215512~231.6631512774GHz"
#cont_channels31 = "31:244~250;323~332;373~377;434~440;467~473;519~533;598~608;623~626;687~695;775~780;1116~1124;1337~1338;1476~1477;1747~1750;1894~1898"
cont_channels31 = "31:232.1391218216~232.1449812745GHz;232.2162712852~232.2250604645GHz;232.2651000596~232.2690063615GHz;232.3246711644~232.3305306173GHz;232.3568981555~232.3627576084GHz;232.4076800809~232.4213521378GHz;232.4848295445~232.4945952994GHz;232.5092439317~232.5121736582GHz;232.571744763~232.5795573669GHz;232.657683406~232.6625662834GHz;232.9906956476~232.9985082515GHz;233.2065188306~233.2074954061GHz;233.3422628235~233.343239399GHz;233.6069147809~233.6098445074GHz;233.7504713778~233.7543776797GHz"

cont_channels = cont_channels25 + ',' + cont_channels27 + ',' + cont_channels29 + ',' + cont_channels31

print("# Subtracting continuum...")
os.system('rm -rf '+ msdata + '.contsub25')
uvcontsub(vis=msdata, outputvis=msdata + '.contsub25', field='"G12.89+0.49"', spw='25', fitspec=cont_channels25, fitorder=1, writemodel=False)
os.system('rm -rf '+ msdata + '.contsub27')
uvcontsub(vis=msdata, outputvis=msdata + '.contsub27', field='"G12.89+0.49"', spw='27', fitspec=cont_channels27, fitorder=1, writemodel=False)
os.system('rm -rf '+ msdata + '.contsub29')
uvcontsub(vis=msdata, outputvis=msdata + '.contsub29', field='"G12.89+0.49"', spw='29', fitspec=cont_channels29, fitorder=1, writemodel=False)
os.system('rm -rf '+ msdata + '.contsub31')
uvcontsub(vis=msdata, outputvis=msdata + '.contsub31', field='"G12.89+0.49"', spw='31', fitspec=cont_channels31, fitorder=1, writemodel=False)

####################################
# Self-calibration 
####################################

print("# Phase self-calibration...")

imname ='targetImageForSelfcal'

make_clean(msdata, 'first_'+imname, mfs=True, inpthreshold='1.8mJy', inpcont_channels=cont_channels)
os.system('rm -rf '+ asdmname + '.phase_int.tb1')
gaincal(vis=msdata,
        caltable = asdmname + '.phase_int.tb1',
        field = '"G12.89+0.49"',
        spw=cont_channels,
        solint='inf',
        calmode='p',
        refant='DV04,DA55,DV10',
        gaintype='G')
applycal(vis=msdata,
         field = '"G12.89+0.49"',
         gaintable=[asdmname + '.phase_int.tb1'],
         flagbackup = False)
os.system('rm -rf '+ msdata + "1*")
split(vis=msdata,
      outputvis=msdata + "1",
      datacolumn='corrected')

make_clean(msdata + "1", 'second_'+imname, mfs=True, inpthreshold='1.3mJy', inpcont_channels=cont_channels)
os.system('rm -rf '+ asdmname + '.phase_int.tb2')
gaincal(vis=msdata + "1",
        caltable = asdmname + '.phase_int.tb2',
        field = '"G12.89+0.49"',
        spw=cont_channels,
        solint='int', # int = 6s
        calmode='p',
        refant='DV04,DA55,DV10',
        gaintype='G')
applycal(vis=msdata + "1",
         field = '"G12.89+0.49"',
         gaintable=[asdmname + '.phase_int.tb2'],
         flagbackup = False)
os.system('rm -rf '+ msdata + "2*")
split(vis=msdata + "1",
      outputvis=msdata + "2",
      datacolumn='corrected')

make_clean(msdata + "2", 'third_'+imname, mfs=True, inpthreshold='1.2mJy', inpcont_channels=cont_channels)
os.system('rm -rf '+ asdmname + '.phase_int.tb3')
gaincal(vis=msdata + "2",
        caltable = asdmname + '.phase_int.tb3',
        field = '"G12.89+0.49"',
        spw=cont_channels,
        solint='int', # int = 6s
        calmode='p',
        refant='DV04,DA55,DV10',
        gaintype='G')
applycal(vis=msdata + "2",
         field = '"G12.89+0.49"',
         gaintable=[asdmname + '.phase_int.tb3'],
         flagbackup = False)
os.system('rm -rf '+ msdata + "3*")
split(vis=msdata + "2",
      outputvis=msdata + "3",
      datacolumn='corrected')

make_clean(msdata + "3", 'fourth_'+imname, mfs=True, inpthreshold='1.2mJy', inpcont_channels=cont_channels)

####################################
# Line imaging
####################################

print("# Line imaging...")

applycal(vis=msdata + '.contsub25',
         field = '"G12.89+0.49"',
         gaintable=[asdmname + '.phase_int.tb3'],
         flagbackup = False)
make_clean(msdata + '.contsub25', 'line25_'+imname, mfs=False, inpthreshold='0.015Jy', inpcont_channels=['25'])

applycal(vis=msdata + '.contsub27',
         field = '"G12.89+0.49"',
         gaintable=[asdmname + '.phase_int.tb3'],
         flagbackup = False)
make_clean(msdata + '.contsub27', 'line27_'+imname, mfs=False, inpthreshold='0.015Jy')

applycal(vis=msdata + '.contsub29',
         field = '"G12.89+0.49"',
         gaintable=[asdmname + '.phase_int.tb3'],
         flagbackup = False)
make_clean(msdata + '.contsub29', 'line29_'+imname, mfs=False, inpthreshold='0.015Jy', inpcont_channels=['29'])

applycal(vis=msdata + '.contsub31',
         field = '"G12.89+0.49"',
         gaintable=[asdmname + '.phase_int.tb3'],
         flagbackup = False)
make_clean(msdata + '.contsub31', 'line31_'+imname, mfs=False, inpthreshold='0.012Jy', inpcont_channels=['31'])
