#!/usr/bin/env python3
import utilities as u
import sys
import os
from subprocess import call
from datetime import datetime,timedelta



#
# for overall limits with a particular resolution - break up into bite size pieces base on maxPoints
#
def findLimits(limits, resolution) :
    sX=int(limits[2]/resolution)
    sY=int(limits[3]/resolution)
    dxkm=resolution
    dykm=resolution
    y0orig=limits[1]
    x0=limits[0]
    gd=u.geodat(x0,y0orig,sX,sY,resolution*1000,resolution*1000)
    # maximum number of points in a file
    maxPoints=8400000
    # this will force more boxes for 0.5 km and greater to run faster
    if resolution > 0.49 :
        maxPoints=maxPoints/4
    # force it to give integer size
    maxLines=int(maxPoints/sX)
    # iteration to adjust size so last case is not too small (>.25)
    done=False
    while not done :
        frac=float(sY)/maxLines - int(sY/maxLines)
        if frac > 0.25 :
            done=True
        else :
            maxLines -=100
    line=0
    lastline=0
    limits=[]
    while line < sY :
        deltaLines=int(min(sY-line,maxLines))
        line += deltaLines
        y0=y0orig+lastline*dykm
        lastline=line
        limits.append([x0,y0,sX*dxkm,deltaLines*dykm,dxkm,dykm])
    return limits,gd
#-----------------------------------------------------------------------------
# Read in input file,  save as list with [date, line] entries. 
#-----------------------------------------------------------------------------
def processInputFile(filename,noTSX) :
    inputFile=open(filename,'r')
    dataTakes=[]
    for line in inputFile :
        # use line length to skip over non data stuff
        if not ';' in line and len(line) > 120 :
            if  noTSX and 'TSX' in line :
                continue
            # setup cross flag, 1 if no TSX, 0 if TSX
            crossFlag=['1','0']['TSX' in line]
            line=line.strip('\n')
            parts=line.split()
            # get geodat
            geo=parts[1].strip()
            exclude="/".join(geo.split('/')[0:-2])+'/Exclude'
            if not os.path.exists(exclude) :
                fgeo=open(geo,'r')
                # read geodat to extract date
                for gline in fgeo :
                    if 'Image date' in gline :
                        gline=gline.split(':')[1].strip('\n').replace(' ','')
                        date=datetime.strptime(gline,'%d%b%Y')
                        break
                fgeo.close()
                line +=' '+crossFlag
                dataTakes.append([date,line])
    inputFile.close()
    return dataTakes

#-----------------------------------------------------------------------------
# Write the new input file, filtering by date.  Pass anything
# that has some overlap with the date range (date1,date2)
# no return value
#-----------------------------------------------------------------------------
def writeInputFile(infile,limit,dataTakes,date1,date2):
    #
    fOut=open(infile,'w')
    print('; automatically generated input file from setupannual mosaics',file=fOut)
    # print limit
    print('{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}'.format(*limit),file=fOut)
    # first pass through to filter by date, and get count
    n=0
    newList=[]
    for dataTake in dataTakes :
        firstDate=dataTake[0]
        tmp=dataTake[1].split()
        secondDate=firstDate+timedelta(float(tmp[3]))
        # pass dates with any amount of overlap
        if  secondDate >= date1 and firstDate <=date2 :
            n+=1
            newList.append(dataTake[1])
    print(len(newList),file=fOut)
    for item in newList :
        print(';\n',item,file=fOut)
#
# Write the geodats for the final mosaic as specified by "geo" in dir vdir
#
def writeGeoDats(vdir,geo) :
    mosaic=vdir+'/'+'mosaicOffsets'
    roots=['.vx.geodat','.vy.geodat','.vz.geodat','.ex.geodat','.ey.geodat']
    for root in roots :
        geoFile=mosaic+root
        geo.writeGeodat(geoFile)

def usage():
    #
    print('\n\t\t\t\t \033[1;34m SETUP ANNUAL VELOCITY MOSAICS \033[0m')
    print('\n'+'\033[1m'+'Use a template inputFile to set up runmosaic.py files to create a series of mosaics, with entries specified')
    print('in a dateLines file (see Notes below). The executable runmosaic.py will take a large mosaic, break it to process it in pieces, runnning ')
    print('multiple threads in parallel. It will then reassemble into a complete mosaic, moving the indivdidual pieces to a subdirectory ./parts\n')
    print("\nsetupannualmosaics.py -dem=demfile -limits='[x0km,y0km,sxkm,xykm]' -baseFlags='-flag1...' -input=inputFile -dates=dateLinesFile -noTSX\n ")
    print('where :\n')
    print('\t -dem\t\t\tis the full path to the DEM file (default: /Volumes/insar3/ian/gimp/gimpdem_90m.270)')
    print('\t -limits\t\tspecifies the output area - use quotes and brackets (default: [-659,-3379,1517,2740])')
    print('\t -baseFlags\t\tmosaic3d flags that are applied to all cases (default: -offsets -rOffsets -fl 50 -noVh)')
    print('\t -inputFile\t\tinputFile with all tracks to be considered (default: track-all/inputFile)')
    print('\t -dates\t\t\tfile with base name, dates, resolution, and flags (default: ./dateLines)')
    print('\t -noTSX\t\t\texclude TSX files (default: include)')    
    print('\nNotes:\n');
    print('Sample dateLines entry - first four required, everything afterwards is a flag :\n')
    print('\tall2016 09-01-2015 05-31-2016 0.5 -shelfMask /Users/ian/greenlandmask/Winter201516/greenlandmask201516 -timeThresh 36 -3dOff\n'+'\033[0m')
    exit()

    
def processArg(argv) :
    dem='/Volumes/insar3/ian/gimp/gimpdem_90m.270'
    limits=[-659,-3379,1517,2740]
    baseFlags='-offsets -rOffsets -fl 50 -noVh'
    inputFile='track-all/inputFile'
    dateLinesFile='dateLines'
    noTSX=False
    try :
        for s in argv[1:] :
            if '-dem' in s :
                dem=s.split('=')[-1]
                print('dem = ',dem)
            elif '-help' in s :
                usage()
            elif '-noTSX' in s :
                noTSX=True                                
            elif '-limits' in s :
                limits=eval(s.split('=')[-1])
                print('limits = ',limits)
            elif '-baseFlags' in s:
                temp=s.split('=')[1:]
                baseFlags=' '.join(temp)
                print('baseFlags = ',baseFlags)
            elif '-input' in s:
                inputFile=s.split('=')[-1]
                print('inputFile = ',inputFile)
            elif '-dates' in s:
                dateLinesFile=s.split('=')[-1]
                print('dataLinesFile = ',dateLinesFile)                
            else :
                usage()
                exit()
        return dem,limits,baseFlags,inputFile,dateLinesFile,noTSX
    except Exception :
        usage()
        raise NameError()
        exit()

     
#
# Set up mosaic files, makes directories if needed, produces all files, can run with runAll
#------------------------------------------------------------------------------------------------------------------------
def main() :
    #
    dem,limits,baseFlags,inputFile,dateLinesFile,noTSX=processArg(sys.argv)
    #
    dataTakes=processInputFile(inputFile,noTSX)
    #
    dateLines=open(dateLinesFile,'r')
    #
    # process each date entry
    #
    for dateLine in dateLines :
        # skip blank lines
        if dateLine.isspace() :
            continue
        # process line
        parts=dateLine.split()
        base=parts[0]
        date1S=parts[1]
        date2S=parts[2]
        resolution=parts[3]
        extraFlags=''
        if len(parts) > 4 :
            extraFlags=' '.join(parts[4:])
        # velocity directiory
        vdir='track-'+base.strip()+'-'+str(resolution)+'km'
        if not os.path.isdir(vdir) :
            os.mkdir(vdir)
        if not os.path.isdir(vdir+'/parts') :
            os.mkdir(vdir+'/parts')
        # open runfile
        runFile=vdir+'/' + 'runmosaic.py'
        fRun=open(runFile,'w')
        print('#!/usr/bin/env python3',file=fRun)
        print('import threading',file=fRun)
        print('import time',file=fRun)
        print('from subprocess import call',file=fRun)
        print('import utilities as u\n#',file=fRun)
        print('dem="',dem,'"','\n#',file=fRun)
        print('threads=[]\n#',file=fRun)
        #
        resolution=float(resolution)
        date1=datetime.strptime(date1S,'%m-%d-%Y')
        date2=datetime.strptime(date2S,'%m-%d-%Y')
        flags=base
        # break it up into to pieces
        pieces,geo=findLimits(limits,resolution)
        #
        writeGeoDats(vdir,geo)
        #
        print(pieces)
        #
        count=0
        # make input files
        for piece in pieces :
            countStr='{:03}'.format(count)
            infile=vdir+'/inputFile-'+countStr
            writeInputFile(infile,piece,dataTakes,date1,date2)
            flags=baseFlags+' '+extraFlags+' -date1 ' + date1.strftime('%m-%d-%Y')+' -date2 ' + date2.strftime('%m-%d-%Y')
            #command='mosaic3d '+flags+ ' inputFile-'+countStr
            print('flags="'+flags+'"',file=fRun)
            print('inputFile="inputFile-'+countStr+'"',file=fRun)
            print('outputFile="mosaicOffsets-'+countStr+'"',file=fRun)
            print('args=[dem,flags,inputFile,outputFile]',file=fRun)
            print('thread'+countStr+'=threading.Thread(target=u.runvel,args=args)',file=fRun)
            print('threads.append(thread'+countStr+')',file=fRun)
            print('#\n# Piece '+countStr+'\n#',file=fRun)
            count+=1
        #
        print('print(threads)',file=fRun)
        #
        # create code to run threads
        #
        print('#\n# Run Threads\n#',file=fRun)
        print('for thread in threads :',file=fRun)
        print('\tthread.start()',file=fRun)
        print('\ttime.sleep(1)',file=fRun)
        print('#\n# Monitor Threads\n#',file=fRun)
        print('nDone=0',file=fRun)
        print('seconds=0',file=fRun)
        print("print('\\n******** RUNNING MOSAIC PIECES *****\\n\\n')"  ,file=fRun)
        print('while nDone < len(threads):',file=fRun)
        print('\tnDone=0',file=fRun)
        print('\tfor thread in threads :',file=fRun)
        print('\t\tif not thread.isAlive() :',file=fRun)
        print('\t\t\tnDone+=1',file=fRun)
        print("\tprint('Time ','{:7.2f}'.format(seconds/60),'(m) N threads completed ',nDone,'N threads running ',len(threads)-nDone,end='\\r')",file=fRun)
        print('\ttime.sleep(5)',file=fRun)
        print('\tseconds+=5',file=fRun)
        print('#\n# Cat file \n#',file=fRun)
        print("call('cat mosaicOffsets-*.vx > mosaicOffsets.vx',shell=True)",file=fRun)
        print("call('cat mosaicOffsets-*.vy > mosaicOffsets.vy',shell=True)",file=fRun)
        print("call('cat mosaicOffsets-*.vz > mosaicOffsets.vz',shell=True)",file=fRun)
        print("call('cat mosaicOffsets-*.ex > mosaicOffsets.ex',shell=True)",file=fRun)
        print("call('cat mosaicOffsets-*.ey > mosaicOffsets.ey',shell=True)",file=fRun)
        print("call('mv mosaicOffsets-* parts',shell=True)",file=fRun)
        print("print('\\n\\n******** END MOSAIC *****\\n')"  ,file=fRun)
        fRun.close()
        call('chmod +x '+runFile,shell=True)
    dateLines.close()
    
main()
