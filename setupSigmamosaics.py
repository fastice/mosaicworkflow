#!/usr/bin/env python3
import utilities as u
import sys
import os
from subprocess import call
from datetime import datetime,timedelta,date
import numpy as np
import stat
import argparse



def setupSigArgs():
    ''' Handle command line args'''
    #print(''\033[1mSetup a CSK HDF for processing\033[0m')
    parser=argparse.ArgumentParser(description='Setup sigma mosaics',
                                   epilog='Notes :  ')
    parser.add_argument('--tracks',nargs='+',default=[26,74,90,112,141,170,83],type=int,
                        help='track numbers [26,74,90,112,141,170,83]')
    parser.add_argument('--tiffOnly',action='store_true',default=False,help='make runAlltiffOnly to just redo tiffs')
    parser.add_argument('--rootDir',nargs=1,default='/Volumes/insar9/ian/Sentinel1',type=str,
                        help='root dir')
    args=parser.parse_args()
    print(args.rootDir,args.tracks)
    return args.rootDir,args.tracks,args.tiffOnly
#--------------------------------------------------------------------------------
#read geodate and return list with info
#--------------------------------------------------------------------------------
def readGeodat(file) :
    geoFile=open(file,'r')
    #
    print(file)
    geoList=[]
    for line in geoFile:
        if not any(s in line for s in [';','#','&']) and not line.isspace() :
            geoList.append(line.split())
    geoData=np.array(geoList).astype(np.float)
    return geoData
#--------------------------------------------------------------------------------
# process geo structure, slice up data, and make headers
#--------------------------------------------------------------------------------
def processGeo(geo) :
    #
    #
    maxBuf=0.3e9
    #
    geodat=geo[1]
    done=False
    #
    # Iterate to check that
    #
    nLines=geodat[0][1]
    nx=geodat[0][0]
    y0=geodat[2][1]
    y0orig=y0
    x0=geodat[2][0]
    dx=geodat[1][0]
    dy=geodat[1][1]
    dxkm=float(geodat[1][0]/1000.)
    dykm=float(geodat[1][1]/1000.)
    #
    # Compute first guess maxlines for buffer
    #
    maxLines=int(maxBuf/(geodat[0][0]*4))
    # force it to give integer size
    maxLines=int(maxLines*(dy/1000))/(dykm)
    # iteration to adjust size so last case is not too small (>.25)
    while not done :
        frac=float(nLines)/maxLines - int(geodat[0][1]/maxLines)
        if frac > 0.25 :
            done=True
        else :
            maxLines -=100
        #
        # now compute limites
        #
        line=0
        lastline=0
        limits=[]
        while line < nLines :
            deltaLines=int(min(nLines-line,maxLines))
            line += deltaLines
            y0=y0orig+lastline*dykm
#            print(x0,y0,nx*dxkm,deltaLines*dykm,dxkm,dykm)
            lastline=line
            limits.append([x0,y0,nx*dxkm,deltaLines*dykm,dxkm,dykm])
        return limits

#--------------------------------------------------------------------------------
# For a given track in rootDir, this routine will :
# a) get a listing of all the geodats - where nl=MxN
# b) extract date info
# c) get pow file
#--------------------------------------------------------------------------------
def getInputsForTrack(track,rootDir,nl) :
    # geodatfiles
    geos=u.dols('ls '+ rootDir+'/track-'+str(track)+'/*_*/geodat'+nl+'.in')
    #
    # Loop through each geo and
    #
    trackData=[]
    for geo in geos :
        rawRoot=geo.split("_")[0]
        # get gain file
        try :
            gainFile=u.dols('ls '+rawRoot+'-*/absolutegain')[0]
        except :
            gainFile=' '
        #
        # NOTE I am assuming gain is same for whole orbit so if multiple sequences this won't be a problem.
        #
        #if os.path.exists(gainFile) :
        #    fpgain=open(gainFile,'r')
        #   line=fpgain.readline()
        #    g1,g2,g3=line.split()
        #    gain=1./((float(g1)+float(g2)+float(g3))/3.)
        #    fpgain.close()
        #else :
        #    print('\n\t\t\033[1;34m  Warning no absolute gain for  \033[0m: ',geo)
        #    gain=1.0
        gain=1.0
        #
        # This block reads the date from the geo file
        #
        fgeo=open(geo,'r')
        for line in fgeo :
            if 'Image date' in line :
                line=line.split(':')[1].strip('\n').replace(' ','')
                date=datetime.strptime(line,'%d%b%Y')
                break
        fgeo.close()
        line1=geo.replace('geodat'+nl+'.in','')
        powFile=u.dols('ls '+line1+'/P*'+nl+'.pow')
        if len(powFile) > 0 :
            powFile=powFile[0]
        # only save files where length > rootDir (assume valid return from dols
        if len(geo) > len(rootDir) and len(powFile) > len(rootDir) :
            result=[date,geo,powFile,gain]
            trackData.append(result)
    return trackData
#--------------------------------------------------------------------------------
# print input files for each track
#-------------------------------------------------------------------------------
def printTracks(allTracks,geoInfo,firstDate,lastDate,midDate,dem,tiffOnly) :
    # direction to make everything in
    dateDir=firstDate.strftime("%Y-%m-%d-sigma0")
    # check
    inputLines=[]
    for track in allTracks :
        if track[0] >= firstDate and track[0] <= lastDate :
            inputLines.append(track)
            validData=True
    # return if no valid data
    if len(inputLines) < 1 :
        print('No Date for '+firstDate.strftime("%Y-%m-%d") + lastDate.strftime("%Y-%m-%d"))
        return False
    #
    # At least some good data, so contine
    #
    runFile=dateDir+'/runAll'
    if tiffOnly :
        runFile +='tiffOnly'
    fRun=open(runFile,'w')
    print('#',file=fRun)
    print('set DEM='+dem,file=fRun)
    #
    print('#\nmakeimageshapefile.py inputFile-NorthEast-0 '+firstDate.strftime('%m-%d-%Y')+'-'+lastDate.strftime('%m-%d-%Y')+'.sigma0',file=fRun)
    print('#\n',file=fRun)
    for region in geoInfo :
        prefix=region[0].split('.geodat')[0]
        print('#\n# ',prefix,'\n#',file=fRun)
        count=0
        filesUsed,filesUsedInc,filesUsedGam='','',''
        for limit in region[2] :
            # add entry to runAll
            if not tiffOnly :
                inputFile='inputFile-'+prefix+'-'+str(count)
                #print(inputFile,limit,end='\n')
                print('# ',count,file=fRun)
                print('geomosaic -center  -S1Cal -smoothOut 2 -descending '+' \\',file=fRun)
                print('\t-date1 '+firstDate.strftime('%m-%d-%Y')+' -date2 '+lastDate.strftime('%m-%d-%Y')+
                      ' -removePad 45 -fl 25 -xyDEM \\',file=fRun)
                print('\t'+inputFile+' $DEM '+prefix+'-'+str(count),file=fRun)
                filesUsed += ' '+prefix+'-'+str(count)
                filesUsedInc += ' '+prefix+'-'+str(count)+'.inc'
                filesUsedGam += ' '+prefix+'-'+str(count)+'.gamcor'
                #
                # now make input file
                #
                finput=open(dateDir+'/'+inputFile,'w')
                # put dates as comment
                print(';\n;',firstDate.strftime('%m-%d-%Y'),' - ',lastDate.strftime('%m-%d-%Y')+'\n;',file=finput)
                # print geo info
                print( *limit ,file=finput)
                # number of items
                print(';',file=finput)
                print(len(inputLines),file=finput)
                print(';',file=finput)
                # output lines
                for line in inputLines :
                    print(line[2],' ',line[1],' ','{:.5}'.format(line[3]),file=finput)
                print('; end',file=finput)
                finput.close()
            count+=1
        # cate files and remove pieces
        outfile=prefix+'.'+firstDate.strftime('%Y-%m-%d')+'.to.'+lastDate.strftime('%Y-%m-%d')
        if not tiffOnly :
            print('#\n# Cat files for '+prefix+'\n#',file=fRun)
            print('cat ' + filesUsed + ' > '+ outfile,file=fRun)
            print('rm ' + filesUsed,file=fRun)
            print('# incidence angle',file=fRun);
            print('cat ' + filesUsedInc + ' > '+ outfile+'.inc',file=fRun)
            print('cp '+outfile+'.geodat ' + outfile+'.inc.geodat',file=fRun)
            print('rm ' + filesUsedInc,file=fRun)
            #
            print('# gamma corr',file=fRun);
            print('cat ' + filesUsedGam + ' > '+ outfile+'.gamcor',file=fRun)
            print('cp '+outfile+'.geodat ' + outfile+'.gamcor.geodat',file=fRun)
            print('rm ' + filesUsedGam,file=fRun)
            #
            print('rm ' + prefix + '-*.geodat',file=fRun)
            # make geodat
            geoFile=dateDir+'/'+outfile+'.geodat'
            fgeo=open(geoFile,'w')
            print('# 2',file=fgeo)
            geo=region[1]
            print(int(geo[0][0]),int(geo[0][1]),file=fgeo)
            print(*geo[1],file=fgeo)
            print(*geo[2],file=fgeo)
            print('&',file=fgeo)
            fgeo.close()
        #idlCommand='idl -e "maketops100tif,'+"'"+outfile+"','"+outfile+".tif'"+'"'
        #print('#\n'+idlCommand,file=fRun)
        #print('# process tile, X X supresses pyramids, and jpgs',file=fRun)
        #print('processTile '+outfile+ ' X X',file=fRun)
        print('#',file=fRun)
        print('sigmatotiff.py --predictor 1 --epsg 3413 --suffix "geo.sigma0" '+outfile,file=fRun)
        print('sigmatotiff.py --predictor 2 --epsg 3413 --suffix ".sigma0" '+outfile+'.inc',file=fRun)
        print('sigmatotiff.py --predictor 1 --epsg 3413 --suffix ".sigma0" '+outfile+'.gamcor',file=fRun)
        #print('rm '+outfile+'.tif',file=fRun)
    print('#\n# Process tiles to create sigma VRT \n#',file=fRun)
    finalMosaic='mosaic.'+ firstDate.strftime('%Y-%m-%d') + '.to.' + lastDate.strftime('%Y-%m-%d')+ '.sigma0.vrt'
    jpgMosaic='mosaic.'+ firstDate.strftime('%Y-%m-%d') + '.to.' + lastDate.strftime('%Y-%m-%d')+ '.sigma0.jpg'
    print('gdalbuildvrt  -srcnodata "-30." -vrtnodata "-30." ' +finalMosaic + ' *geo.sigma0.tif',file=fRun)
    print('gdaladdo -clean ' +finalMosaic,file=fRun)
    print('gdaladdo --config COMPRESS_OVERVIEW LZW --config TILED YES -r average ' +finalMosaic+' 32 16 8 4 2',file=fRun)
    print('gdalinfo -stats ' +finalMosaic,file=fRun)
    print('gdal_translate -outsize 5% 5% -ot Byte -scale -30 20  -of jpeg  '+finalMosaic+' ' + jpgMosaic ,file=fRun)
    # incidence angle
    print('# Incidence Angle',file=fRun);
    finalMosaicInc='mosaic.'+ firstDate.strftime('%Y-%m-%d') + '.to.' + lastDate.strftime('%Y-%m-%d')+ '.inc.vrt'
    print('gdalbuildvrt  -srcnodata "0." -vrtnodata "0." ' +finalMosaicInc + ' *.inc.sigma0.tif',file=fRun)
    print('gdaladdo --config COMPRESS_OVERVIEW LZW --config TILED YES -r average ' +finalMosaicInc +' 32 16 8 4 2',file=fRun)
    print('gdalinfo -stats ' +finalMosaicInc,file=fRun)
    #
    print('# gamma corr',file=fRun);
    finalMosaicGam='mosaic.'+ firstDate.strftime('%Y-%m-%d') + '.to.' + lastDate.strftime('%Y-%m-%d')+ '.gamcor.vrt'
    print('gdalbuildvrt  -srcnodata "-30." -vrtnodata "-30." ' +finalMosaicGam + ' *.gamcor.sigma0.tif',file=fRun)
    print('gdaladdo --config COMPRESS_OVERVIEW LZW --config TILED YES -r average ' +finalMosaicGam +' 32 16 8 4 2',file=fRun)
    print('gdalinfo -stats ' +finalMosaicGam,file=fRun)
    print('#----------------------------------------done-------------------------------------------------------\n',file=fRun)
    fRun.close()
    os.chmod(runFile,os.stat(runFile).st_mode | stat.S_IEXEC)
    return True

#
# Set up mosaic files, makes directories if needed, produces all files, can run with runAll
#------------------------------------------------------------------------------------------------------------------------
def main() :
    #dem='/Volumes/insar7/ian/fromInsar3/gimp/version3/50m/gimpversion3.0.50'
    dem='/Volumes/insar7/ian/gimp/gimp2/50m/dem.gimp2.50m'
    nl='10x2'
    #
    # process geodat inf an seg up with limtes
    #
    rootDir,tracks,tiffOnly=setupSigArgs()
    u.myerror('Obsolete - use makeimagemosaics.py')
    geodats=u.dols('ls *.geodat')
    if len(geodats) < 4 :
        print('Error should have NorthEast.geodat, North.geodat, NorthWest.geodat, South.geodat as templates')
    geoInfo=[]
    for geodat in geodats :
        geoInfo.append([geodat,readGeodat(geodat)])
    # stores all geoinfo in lists [file,geodat data, list of resolution data]
    for geo in geoInfo :
        limits=processGeo(geo)
        geo.append(limits)
    #

    dateLines=open('sendates','r')
    #
    # getinput info
    #

    allTracks=[]
    for track in tracks :
        print('Processing track '+str(track))
        try :
            oneTrackDat=getInputsForTrack(track,rootDir,nl)
            allTracks += oneTrackDat
        except NameError :
            print('Problem parsing track '+str(track))
    #
    # loop through each date
    #
    for dateLine in dateLines :
        # skip this iteration if no data
        if dateLine.isspace() :
            continue
        #
        # extract dates as python dates
        #
        firstDate,lastDate,midDate = dateLine.strip('\n').split()
        firstDate=datetime.strptime(firstDate,'%m-%d-%Y')
        lastDate=datetime.strptime(lastDate,'%m-%d-%Y')
        midDate=datetime.strptime(midDate,'%m-%d-%Y')
        #
        # form date dir name, check existence, make dir if none
        #
        dateDir=firstDate.strftime('%Y-%m-%d-sigma0')
        createdDir=False
        if not os.path.isdir(dateDir) :
            print('creating '+dateDir)
            os.mkdir(dateDir)
            createdDir=True
        #
        # cycle through geodats
        #
        tracksPrinted=printTracks(allTracks,geoInfo,firstDate,lastDate,midDate,dem,tiffOnly)
        # if just made directory, but nothing printed, remove the new directory
        if not tracksPrinted and createdDir :
            os.rmdir(dateDir)

    dateLines.close()

main()
